#!/usr/bin/python 
import os, sys, threading
from scipy.optimize import minimize

CUBE_CHARGE_SCALE = 100

class Embedding:
	def __init__(self):
		self.symbols = []
		self.limits = []
		self.groups = {}
		self.cube = []

	def full_value(self, ls):
		res = []
		n = 0
		for i, s in enumerate(self.symbols):
			found = False
			for k in self.groups.keys():
				if i in self.groups[k][1:]:
					root = self.groups[k][0]
					res.append(res[root])
					found = True
			if not found:
				res.append(ls[n])
				n += 1
		return res

	def group_list(self):
		res = []
		for i, s in enumerate(self.symbols):
			if i in self.cube:
				res.append("CUBE")
			else:
				found = False
				for k in self.groups.keys():
					if i in self.groups[k]:
						res.append(k)
						found = True
				if not found:
					res.append("")
		return res

	def full_charges(self, x):
		chgs = self.full_value(x)
		for i, c in enumerate(chgs):
			if i in self.cube:
				chgs[i] *= CUBE_CHARGE_SCALE
		return chgs

	def cube_sum(self, x):
		chgs = self.full_charges(x)
		res = 0
		for i, c in enumerate(chgs):
			if i in self.cube:
				res += c
		return res

	def full_limits(self):
		lmts = self.full_value(self.limits)
		for i, (l, h) in enumerate(lmts):
			if i in self.cube:
				if h != None:
					h *= CUBE_CHARGE_SCALE
				if l != None:
					l *= CUBE_CHARGE_SCALE
				lmts[i] = (l, h)
		return lmts

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
			e.symbols.append(ls[0])
			float_or_none = lambda x: float(x) if x != '*' else None
			if len(ls) == 4:
				charges.append(float(ls[1]))
				e.limits.append((float_or_none(ls[2]), float_or_none(ls[3])))
			elif len(ls) == 5:
				charge = float(ls[1])
				limits = (float_or_none(ls[2]), float_or_none(ls[3]))
				grname = ls[4]
				if grname == 'CUBE':
					charges.append(charge / CUBE_CHARGE_SCALE)
					e.limits.append([x / CUBE_CHARGE_SCALE if x != None else x for x in limits])
					e.cube.append(i)
				elif grname not in e.groups.keys():
					charges.append(charge)
					e.limits.append(limits)
					e.groups[grname] = [i]
				else:
					root = e.groups[grname][0]
					assert charge == charges[root], "All charges for same group must be equal"
					assert limits == e.limits[root], "All limits for same group must be equal"
					e.groups[grname].append(i)
			else:
				raise Exception("That's a very bad line: %s" % l)
	return e, charges

def read_embed_positions():
	res = []
	with open("coord") as c:
		n = 0
		for l in c:
			ls = l.split()
			if '$' not in l and ls:
				if ls[-1] == 'zz':
					res += [n]
				n += 1
	return n, res

def read_grad_from_control(num, pos):
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
				if i not in pos:
					rawgrads.extend(igrads)
			grad = (reduce(lambda s, i: s + i**2, rawgrads, 0.))**0.5

def write_embedding(e, x):
	with open("embedding", "w") as emb:
		write_embedding_block(emb, e, x, False)

def write_restart(e, x):
	with open("embedding.restart", "w") as emb:
		write_embedding_block(emb, e, x, True)

def write_result(e, x, grad):
	with open("embedding.log", "a") as emb:
		emb.write("********** GRAD = %13.10f ********\n" % (grad))
		write_embedding_block(emb, e, x, False)

def write_embedding_block(emb, e, x, include_limits_and_groups):
	limit_str = lambda x: "%7.3f" % x if x != None else '*'
	for s, c, (l, h), g in zip(e.symbols, e.full_charges(x), e.full_limits(), e.group_list()):
		if include_limits_and_groups:
			emb.write("%2s %12.8f %7s %7s %s\n" % (s, c, limit_str(l), limit_str(h), g))
		else:
			emb.write("%2s %12.8f\n" % (s, c))

lock = threading.Lock()

def calculate(n, pos, e, x):
	with lock:
		write_embedding(e, x)
		assert os.system('dscf > log_dscf.log') == 0
		assert os.system('grad > log_grad.log') == 0
		grad = read_grad_from_control(n, pos)
		write_result(e, x, grad)
		return grad

if len(sys.argv) == 3:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])

	n, pos = read_embed_positions()
	e, x0 = read_start_embedding()
	s = sum(e.full_charges(x0))
	f = lambda x: calculate(n, pos, e, x)
	cons = ({'type': 'eq', 'fun' : lambda x: sum(e.full_charges(x)) - s}, {'type': 'eq', 'fun' : lambda x: e.cube_sum(x)})

	with open("embedding.log", "w") as emb:
		emb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))

	write_restart(e, x0)
	res = minimize(f, x0, args=(), method='SLSQP', jac=None, 
bounds=e.limits, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': 9001, 'ftol': ftol})
	write_embedding(e, res.x)
	write_restart(e, res.x)
else:
	print 'Usage: charge_optimizer.py [STEP] [PRECISION]\nExample: charge_optimizer.py 1e-2 1e-5'




























