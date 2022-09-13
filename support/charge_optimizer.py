#!/usr/bin/python 
import os, sys, threading, time
from scipy.optimize import minimize
import numpy as np
from functools import reduce


last_xe = None
last_xe_step = None

PREC = 2e-8

def check_if_step_not_grad(xe, eps):
	xe = list(xe)
	global last_xe_step
	global last_xe
	if last_xe == None:
		last_xe_step = xe
		last_xe = xe
		return True
	else:
		assert len(xe) == len(last_xe_step)
		n_diff = 0
		s_diff = 0
		for x, lx in zip(xe, last_xe_step):
			if abs(x - lx) >= PREC:
				n_diff += 1
				s_diff += abs(x - lx)
				
		if n_diff == 1 and abs(s_diff - eps) <= PREC:
			last_xe = xe
			return False
		else:
			last_xe_step = xe
			last_xe = xe
			return True

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
		print("restarting")
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

def read_grad_from_control(num, empos):
	with open("control") as c:
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

def write_embedding(e, xe):
	with open("embedding", "w") as emb:
		write_embedding_block(emb, e, xe, False)

def write_restart(e, xe):
	with open("embedding.restart", "w") as emb:
		write_embedding_block(emb, e, xe, True)

def write_result(e, xe, grad, truestep):
	with open("embedding.log", "a") as emb:
		emb.write("********** VALUE = %12.10f | SUM = %g | TYPE = %s **********\n" % (grad, sum(e.apply_groups(xe)), "MAIN" if truestep else "GRAD"))
		if truestep:
			write_embedding_block(emb, e, xe, False)

def write_embedding_block(emb, e, xe, include_limits_and_groups):
	limit_str = lambda x: "%9.5f" % x if x != None else '*'
	charges = e.apply_groups(xe)
	hlimits = e.apply_groups(e.highlimits)
	llimits = e.apply_groups(e.lowlimits)
	for s, c, l, h, g in zip(e.symbols, charges, llimits, hlimits, e.group_list()):
		if include_limits_and_groups:
			emb.write("%2s %12.8f %10s %10s %s\n" % (s, c, limit_str(l), limit_str(h), g))
		else:
			emb.write("%2s %12.8f\n" % (s, c))

def calculate(n, empos, e, xe, eps):
		write_embedding(e, xe)
		assert os.system('/home/demidov/turbo/bin/amd64/part_dscf > log_dscf.log') == 0
		time.sleep(5)
		assert os.system('/home/demidov/turbo/bin/amd64/part_grad > log_grad.log') == 0
		time.sleep(5)
		grad = read_grad_from_control(n, empos)
		t = check_if_step_not_grad(xe, eps)
		write_result(e, xe, grad, t)
		if t:
			write_restart(e, xe)
		return grad


def optimize_embedding(eps, ftol, maxit):
	e = read_start_embedding()
	x0 = e.charges
	n, empos = read_embed_positions()
	s = sum(e.apply_groups(x0))
	cons = ({'type': 'eq', 'fun' : lambda x: sum(e.apply_groups(x)) - s})
	f = lambda x: calculate(n, empos, e, x, eps)
	limits = list(zip(e.lowlimits, e.highlimits))
	res = minimize(f, x0, args=(), method='SLSQP', jac=None, 
bounds=limits, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
	write_embedding(e, res.x)
	write_restart(e, res.x)
	return res.status, res.nit

if len(sys.argv) == 3:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])
	with open("embedding.log", "w") as emb:
		emb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))
	optimize_embedding(eps, ftol, 9999)

else:
	print('Usage: charge_optimizer.py [STEP] [PRECISION]\nExample: charge_optimizer.py 1e-2 1e-5')




























