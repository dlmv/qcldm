#!/usr/bin/python 
import os, sys, threading
from scipy.optimize import minimize

#CUBE_CHARGE_SCALE = 500

class Embedding:
	def __init__(self):
		self.symbols = []
		self.limits = []
		self.groups = {}
		self.cube = []

	def full_value(self, ls_emb, ls_cub):
		res = []
		ne = 0
		nc = 0
		for i, s in enumerate(self.symbols):
			if self.cube[-1] == i:
				res.append(None)
			elif i in self.cube:
				res.append(ls_cub[nc])
				nc += 1
			else:
				found = False
				for k in self.groups.keys():
					if i in self.groups[k][1:]:
						root = self.groups[k][0]
						res.append(res[root])
						found = True
				if not found:
					res.append(ls_emb[ne])
					ne += 1
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

	def full_charges(self, xe, xc):
		chgs = self.full_value(xe, xc)
		chgs[self.cube[-1]] = -sum(xc)
		return chgs

	def full_limits(self):
		cl = [(None, None)] * (len(self.cube))
		lmts = self.full_value(self.limits, cl)
		lmts[self.cube[-1]] = (None, None)
		return lmts

def read_start_embedding():
	e = Embedding()
	charges = []
	cube_charges = []
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
					cube_charges.append(charge)
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
	if len(cube_charges) > 1:
		cube_charges = cube_charges[:-1]
	return e, charges, cube_charges

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

def write_embedding(e, xe, xc):
	with open("embedding", "w") as emb:
		write_embedding_block(emb, e, xe, xc, False)

def write_restart(e, xe, xc):
	with open("embedding.restart", "w") as emb:
		write_embedding_block(emb, e, xe, xc, True)

def write_result(e, xe, xc, grad):
	with open("embedding.log", "a") as emb:
		emb.write("********** GRAD = %12.10f | SUM = %g **********\n" % (grad, sum(e.full_charges(xe, xc))))
		write_embedding_block(emb, e, xe, xc, False)

def write_embedding_block(emb, e, xe, xc, include_limits_and_groups):
	limit_str = lambda x: "%7.3f" % x if x != None else '*'
	for s, c, (l, h), g in zip(e.symbols, e.full_charges(xe, xc), e.full_limits(), e.group_list()):
		if include_limits_and_groups:
			emb.write("%2s %12.8f %7s %7s %s\n" % (s, c, limit_str(l), limit_str(h), g))
		else:
			emb.write("%2s %12.8f\n" % (s, c))

lock = threading.Lock()

def calculate(n, pos, e, xe, xc):
	with lock:
		write_embedding(e, xe, xc)
		assert os.system('dscf > log_dscf.log') == 0
		assert os.system('grad > log_grad.log') == 0
		grad = read_grad_from_control(n, pos)
#		grad = 0
		write_result(e, xe, xc, grad)
		return grad


def optimize_embedding(e, pos, xe0, xc0, eps, ftol, maxit):
	s = sum(e.full_charges(xe0, xc0))
	f = lambda x: calculate(n, pos, e, x, xc0)
	cons = ({'type': 'eq', 'fun' : lambda x: sum(e.full_charges(x, xc0)) - s})
	res = minimize(f, xe0, args=(), method='SLSQP', jac=None, 
bounds=e.limits, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
	write_embedding(e, res.x, xc0)
	write_restart(e, res.x, xc0)
	return res.status, res.nit

def optimize_cube(e, pos, xe0, xc0, ftol, maxit):
	s = sum(e.full_charges(xe0, xc0))
	f = lambda x: calculate(n, pos, e, xe0, x)
	cons = ()
	res = minimize(f, xc0, args=(), method='SLSQP', jac=None, 
bounds=None, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': 5, 'maxiter': maxit, 'ftol': ftol})
	write_embedding(e, xe0, res.x)
	write_restart(e, xe0, res.x)
	return res.status, res.nit

if len(sys.argv) == 3:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])
	with open("embedding.log", "w") as emb:
		emb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))
	n, pos = read_embed_positions()
	e, xe0, xc0 = read_start_embedding()
	while True:
		st, it = optimize_embedding(e, pos, xe0, xc0, eps, ftol, 5)
		st, it = optimize_cube(e, pos, xe0, xc0, ftol, 5)
		if st == 0 and it == 1:
			break

else:
	print 'Usage: charge_optimizer.py [STEP] [PRECISION]\nExample: charge_optimizer.py 1e-2 1e-5'




























