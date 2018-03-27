#!/usr/bin/python 
import re, os, sys, threading, math, time
from scipy.optimize import minimize

SHELLS = "spdfgh"

class Exponent:
	def __init__(self, c, n, a):
		self.c = c
		self.n = n
		self.a = a

class ECP_Component:
	def __init__(self):
		self.name = ''
		self.exps = []

	def block(self):
		res = self.name + '\n'
		for e in self.exps:
			res += "%20.7f%5d%20.7f\n" % (e.c, e.n, e.a)
		return res
		
		

class Ecp:
	def __init__(self):
		self.name = ''
		self.ndict = {}
		self.components = []
		self.so_components = []

	def block(self):
		res = '*\n'
		res += self.name + '\n*\n'
		for k in ['ncore', 'lmax', 'lsomax']:
			if k in self.ndict.keys():
				res += '  %s=  %d  ' % (k, self.ndict[k])
		res += '\n#        coefficient   r^n          exponent\n'
		for c in self.components:
			res += c.block()
		if self.so_components:
			res += '*\n'
			for c in self.so_components:
				res += c.block()
		return res

class Basis_Data():
	def __init__(self):
		self.basis = ''
		self.ecps = []

def read_component(lines, n):
	res = ECP_Component()
	res.name = lines[n]
	n += 1
	while re.match('\s+.*', lines[n]):
		ls = lines[n].split()
		exp = Exponent(float(ls[0]), int(ls[1]), float(ls[2]))
		res.exps.append(exp)
		n += 1
	return n, res
	
def read_ecp(lines, n):
	res = Ecp()
	n += 1
	res.name = lines[n]
	n += 2
	ls = re.split("[\s=]+", lines[n].strip())
	for k, v in zip(ls[0::2], ls[1::2]):
		res.ndict[k] = int(v)
	n += 2
	while lines[n] != '*':
		n, comp = read_component(lines, n)
		res.components.append(comp)
	if 'lsomax' in res.ndict.keys() and res.ndict['lsomax'] > 0:
		n += 1
		while lines[n] != '*':
			n, comp = read_component(lines, n)
			res.so_components.append(comp)
	return n, res
	
	
def read_basis_and_ecp():
	with open('basis') as bs:
		BD = Basis_Data()
		lines = bs.read().splitlines()
		basis = ''
		n = 0
		while  n < len(lines) and lines[n] != '$ecp':
			basis += lines[n] + '\n'
			n += 1
		BD.basis = basis
		n += 1
		while lines[n + 1] != '$end':
			n, ecp = read_ecp(lines, n)
			BD.ecps.append(ecp)
		return BD

def write_basis_and_ecp(BD):
	with open('basis', 'w') as bs:
		bs.write(BD.basis)
		bs.write('$ecp\n')
		for ecp in BD.ecps:
			bs.write(ecp.block())
		bs.write('*\n$end\n')

def load_var_exponents(ecps):
	res = []
	for ecp in ecps:
		if ecp.name.split()[1].startswith('var_'):
			for comp in ecp.components:
				for exp in comp.exps:
					res.append(math.log10(exp.a))
	return res

def rewrite_var_exponents(ecps, x):
	n = 0
	for ecp in ecps:
		if ecp.name.split()[1].startswith('var_'):
			for comp in ecp.components:
				for exp in comp.exps:
					exp.a = math.pow(10, x[n])
					n += 1
def read_embed_positions():
	pos = []
	nepos = []
	with open("coord") as c:
		n = 0
		for l in c:
			ls = l.split()
			if '$' not in l and ls:
				if ls[-1] == 'zz':
					pos += [n]
				elif ls[-1] == 'ne':#FIXME
					nepos += [n]
				n += 1
	return n, pos, nepos

def read_grad_from_control(num, pos, nepos):
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
				if i not in pos + nepos:
					rawgrads.extend(igrads)
			grad = (reduce(lambda s, i: s + i**2, rawgrads, 0.))**0.5

def write_result(x, grad):
	with open("ecp.log", "a") as emb:
		emb.write("********** GRAD = %12.10f **********\n" % (grad))
		for xx in x:
			emb.write('%20.7f\n' % math.pow(10, xx))


lock = threading.Lock()

def calculate(n, pos, nepos, BD, x):
	with lock:
		rewrite_var_exponents(BD.ecps, x)
		write_basis_and_ecp(BD)
		assert os.system('dscf > log_dscf.log') == 0
		time.sleep(5)
		assert os.system('grad > log_grad.log') == 0
		time.sleep(5)
		grad = read_grad_from_control(n, pos, nepos)
		write_result(x, grad)
		return grad


def optimize_ecp(eps, ftol, maxit):
	BD = read_basis_and_ecp()
	x0 = load_var_exponents(BD.ecps)
	write_basis_and_ecp(BD)
	n, pos, nepos = read_embed_positions()
	f = lambda x: calculate(n, pos, nepos, BD, x)
	res = minimize(f, x0, args=(), method='SLSQP', jac=None, bounds=(), constraints=(), tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
	rewrite_var_exponents(BD.ecps, res.x)
	write_basis_and_ecp(BD)

if len(sys.argv) == 3:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])
	with open("ecp.log", "w") as emb:
		emb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))
	optimize_ecp(eps, ftol, 9999)
else:
	print 'Usage: charge_optimizer.py [STEP] [PRECISION]\nExample: charge_optimizer.py 1e-2 1e-5'




























