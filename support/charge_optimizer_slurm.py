#!/usr/bin/python 
import os, sys, re, time, subprocess, shutil
from scipy.optimize import minimize
import numpy as np



dscf_command = '%s/bin/amd64/part_dscf>log_dscf.log' % os.environ['TURBO']
grad_command = '%s/bin/amd64/part_grad>log_grad.log' % os.environ['TURBO']
srun_command = 'srun -N 1 -n 1 -c 28 %s'

class Embedding:
	def __init__(self):
		self.symbols = []
		self.lowlimits = []
		self.highlimits = []
		self.group_matrix = []
		self.group_names = []
		self.charges = []

	def expand_matrix(self):
		for row in self.group_matrix:
			row.append(0)

	def add_new(self, charge, ll, hl, group):
		self.charges.append(charge)
		self.lowlimits.append(ll)
		self.highlimits.append(hl)
		self.group_names.append(group)
		row = [0] * (len(self.group_matrix[0]) if self.group_matrix else 1)
		row [-1] = 1
		self.group_matrix.append(row)

	def add(self, symbol, charge, ll, hl, group):
		self.symbols.append(symbol)
		self.expand_matrix()
		if group:
			if group in self.group_names:
				index = self.group_names.index(group)
				assert hl == self.highlimits[index], 'upper limits must be equal for equal atoms\nand they are not for group %s:\n%f != %f' % (group, hl, self.highlimits[index])
				assert ll == self.lowlimits[index], 'limits must be equal for equal atoms\nand they are not for group %s:\n%f != %f' % (group, ll, self.lowlimits[index])
				assert charge == self.charges[index], 'charges must be equal for equal atoms\nand they are not for group %s:\n%f != %f' % (group, charge, self.charges[index])
				self.group_matrix[index][-1] = 1
			else:
				self.add_new(charge, ll, hl, group)
		else:
			self.add_new(charge, ll, hl, '')

	def apply_groups(self, x):
		assert len(x) == len(self.charges)
		a = np.array(self.group_matrix).transpose()
		b = np.array(x)
		return list(a.dot(b))

	def group_list(self):
		res = []
		a = np.array(self.group_matrix).transpose()
		for (x,y), value in np.ndenumerate(a):
			if value == 1:
				res.append(self.group_names[y])
		return res
			

def float_or_none(x):
	return float(x) if x != '*' else None

def read_start_embedding():
	e = Embedding()
	charges = []
	filename = "embedding.start"
	if os.path.exists("embedding.restart"):
		filename = "embedding.restart"
		print "restarting"
	with open(filename) as es:
		for i, l in enumerate(es):
			ls = l.split()
			symbol = ls[0]
			charge = float(ls[1])
			ll = float_or_none(ls[2])
			hl =float_or_none(ls[3])
			group = ''
			if len(ls) == 5:
				group = ls[4]
			e.add(symbol, charge, ll, hl, group)	
		return e

def read_embed_positions():
	empos = []
	with open("coord") as c:
		n = 0
		for l in c:
			ls = l.split()
			if '$' not in l and ls:
				if ls[-1] == 'zz':
					empos += [n]
				n += 1
	return n, empos	

def read_grad_from_control(num, empos, folder='.'):
	with open(os.path.join(folder, "control")) as c:
		lines = c.read().splitlines()
		n = 0
		while n < len(lines):
			if '$grad' in lines[n]:
				break
			n += 1
		assert n < len(lines)
		grad = 0.
		while n < len(lines):
			while n < len(lines):
				if lines[n][0] == '#':
					n += 1
					continue
				if 'cycle' in lines[n]:
					break
				n += 1
			if n == len(lines):
				return grad
			grad = 0.
			rawgrads = []
			n0 = n
			for i in range(num):
				n = n0 + num + 1 +  i
				igrads = [float(l.replace("D", "E")) for l in lines[n].split()]
				if i not in empos:
					rawgrads.extend(igrads)
			grad = (reduce(lambda s, i: s + i**2, rawgrads, 0.))**0.5

def write_embedding(e, x, folder='.'):
	with open(os.path.join(folder, "embedding"), "w") as emb:
		write_embedding_block(emb, e, x, False)

def write_restart(e, x):
	with open("embedding.restart", "w") as emb:
		write_embedding_block(emb, e, x, True)

def write_result(e, x, grad, truestep):
	with open("embedding.log", "a") as emb:
		emb.write("********** VALUE = %12.10f | SUM = %g | TYPE = %s **********\n" % (grad, sum(e.apply_groups(x)), "MAIN" if truestep else "GRAD"))
		if truestep:
			write_embedding_block(emb, e, x, False)

def write_embedding_block(emb, e, x, include_limits_and_groups):
	limit_str = lambda x: "%9.5f" % x if x != None else '*'
	charges = e.apply_groups(x)
	hlimits = e.apply_groups(e.highlimits)
	llimits = e.apply_groups(e.lowlimits)
	for s, c, l, h, g in zip(e.symbols, charges, llimits, hlimits, e.group_list()):
		if include_limits_and_groups:
			emb.write("%2s %12.8f %10s %10s %s\n" % (s, c, limit_str(l), limit_str(h), g))
		else:
			emb.write("%2s %12.8f\n" % (s, c))

def all_grads(xs, eps):
	for i, x in enumerate(xs):
		gradlist = list(xs)
		gradlist[i] = x + eps
		yield gradlist
		
class GradTask:
	calc_files = ['coord', 'control', 'realmos', 'imagmos', 'basis']
	def __init__(self, dirname, e, x, cache, gradfunc):
		self.dirname = dirname
		self.e = e
		self.x = x
		self.process = None
		self.state = 0
		self.cache = cache
		self.gradfunc = gradfunc
		
	def create_input(self):
		commands = []
		if os.path.exists(self.dirname):
			shutil.rmtree(self.dirname)
		os.mkdir(self.dirname)
		write_embedding(self.e, self.x, self.dirname)
		for cf in self.calc_files:
			tmpc = "cp %s %s" % (cf, self.dirname)
			commands.append(srun_command % tmpc)
		command = " && ".join(commands)
		self.process = subprocess.Popen(command, shell=True, cwd='.')
		self.state = 1

	def create_input_sync(self):
		commands = []
		if os.path.exists(self.dirname):
			shutil.rmtree(self.dirname)
		os.mkdir(self.dirname)
		write_embedding(self.e, self.x, self.dirname)
		for cf in self.calc_files:
			tmpc = "cp %s %s" % (cf, self.dirname)
			commands.append(tmpc)
		command = " && ".join(commands)
		os.system(command)
		self.state = 1

	def start_dscf(self):
		command = srun_command % dscf_command
		self.process = subprocess.Popen(command, shell=True, cwd=self.dirname)
		self.state = 2

	def start_grad(self):
		command = srun_command % grad_command
		self.process = subprocess.Popen(command, shell=True, cwd=self.dirname)
		self.state = 3

	def status(self):
		if self.process == None:
			return None
		else:
			poll = self.process.poll()
			if poll == None:#working
				return poll
			elif poll != 0:#error
				return poll
			elif self.state == 1:#directory_created:
				self.start_dscf()
				return None
			elif self.state == 2:#dscf finished, but not grad:
				self.start_grad()
				return None
			else:#finished
				value = self.gradfunc(self.dirname)
				self.cache[tuple(list(self.x))] = value
				shutil.rmtree(self.dirname)
				return poll
				

class GradQueue:
	def __init__(self, tasks, limit):
		self.waitq = list(tasks)
		print 'make dirs for all tasks'
		for task in tasks:
			task.create_input_sync()
		self.runq = []
		self.limit = limit

	def try_start_tasks(self):
		while len(self.runq) < self.limit and len(self.waitq) > 0:
			task = self.waitq.pop()
			self.runq.append(task)
			task.start_dscf()

	def check_finished(self, task):
		status = task.status()
		if status is None:
			return False
		elif status == 0:
			return True
		assert False, str(status)

	def update_tasks(self):
#		print 'updating tasks: %d waiting; %d working' % (len(self.waitq), len(self.runq))
		self.runq[:] = [x for x in self.runq if not self.check_finished(x)]
		self.try_start_tasks()

	def do_tasks(self):
		while self.waitq or self.runq:
			self.update_tasks()
			time.sleep(5)

class GradManager:

	def __init__(self, e, n, empos, eps, limit):
		self.cache = {}
		self.last_x = None
		self.last_main_x = None
		self.PREC = 2e-8
		self.grad_started = False
		self.grad_finished = False
		self.e = e
		self.empos = empos
		self.n = n
		self.eps = eps
		self.limit = limit
		
		self.grad_done = False
	
	def check_diff(self, x1, x2):
		assert len(x1) == len(x2)
		n_diff = 0
		s_diff = 0
		for xx1, xx2 in zip(x1, x2):
			if abs(xx1 - xx2) >= self.PREC:
				n_diff += 1
				s_diff += abs(xx1 - xx2)
		return n_diff, s_diff

	def find_cached_value(self, x):
		for k in self.cache.keys():
			n_diff, s_diff = self.check_diff(x, k)
			if n_diff == 0:
				return self.cache[k]
		return None

	def calculate(self, x):
			t = self.check_if_step_not_grad(x, self.eps)
			if t:
				self.grad_done = False
				v = self.find_cached_value(x)
				if v is None:
					v = self.calculate_direct(x)
					self.cache[tuple(x)] = v
				write_result(self.e, x, v, t)
				if t:
					write_restart(self.e, x)
				return v
			else:
				if self.grad_done:
					v = self.find_cached_value(x)
					assert v
					write_result(self.e, x, v, t)
					return v
				else:
					tasks = []
					gradfunc = lambda x: read_grad_from_control(self.n, self.empos, x)
					for i, xs in enumerate(all_grads(self.last_main_x, eps)):
						dirname = 'grad%s' % i
						task = GradTask(dirname, self.e, xs, self.cache, gradfunc)
						tasks.append(task)
					q = GradQueue(tasks, self.limit)
					q.do_tasks()
					self.grad_done = True
					v = self.find_cached_value(x)
					assert v
					write_result(self.e, x, v, t)
					return v
						
			
	
	def calculate_direct(self, x):
		write_embedding(self.e, x)
		assert os.system(srun_command % dscf_command) == 0
		time.sleep(1)
		assert os.system(srun_command % grad_command) == 0
		time.sleep(1)
		grad = read_grad_from_control(self.n, self.empos)
		return grad

	def check_if_step_not_grad(self, x, eps):
		x = list(x)
		if self.last_x == None:
			self.last_main_x = x
			self.last_x = x
			return True
		else:
			n_diff, s_diff = self.check_diff(x, self.last_main_x)
			if n_diff == 1 and abs(s_diff - eps) <= self.PREC:
				self.last_x = x
				return False
			else:
				self.last_main_x = x
				self.last_x = x
				return True

		
def optimize_embedding(eps, ftol, limit, maxit):
	e = read_start_embedding()
	x0 = e.charges
	n, empos = read_embed_positions()
	s = sum(e.apply_groups(x0))
	cons = ({'type': 'eq', 'fun' : lambda x: sum(e.apply_groups(x)) - s})
	G = GradManager(e, n, empos, eps, limit)
	f = lambda x: G.calculate(x)
	limits = zip(e.lowlimits, e.highlimits)
	res = minimize(f, x0, args=(), method='SLSQP', jac=None, 
bounds=limits, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
	write_embedding(e, res.x)
	write_restart(e, res.x)
	return res.status, res.nit

if len(sys.argv) == 4:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])
	limit = int(sys.argv[3])
	with open("embedding.log", "w") as emb:
		emb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))
	optimize_embedding(eps, ftol, limit, 9999)

else:
	print 'Usage: charge_optimizer_slurm.py [STEP] [PRECISION] [LIMIT]\nExample: charge_optimizer.py 1e-2 1e-5 10'




	










