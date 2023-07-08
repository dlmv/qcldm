import re, sys, math
from functools import reduce
sys.dont_write_bytecode = True

from ..numeric_functions.numeric_function import NumericFunction, NumericOperations

def dfac(n):
	return 1 if n < 2 else reduce(lambda x,y: y*x, list(range(n,1,-2)))

def gauss_norm(a, l):
	return (2**(2*l+3.5)  / dfac(2*l+1) / math.pi**0.5)**0.5 * a**((2.*l+3)/4)


class GaussFunction:
	def __init__(self, a, l):
		self.a = a
		self.l = l

	def __call__(self, r):
		return (r**self.l) * math.exp(-self.a * r**2)

class GaussFunctionNormed(GaussFunction):
	def __init__(self, a, l):
		GaussFunction.__init__(self, a,l)
		self.norm = gauss_norm(a, l)

	def __call__(self, r):
		return GaussFunction.__call__(self, r) * self.norm

	@staticmethod
	def overlap(f1, f2):
		if f1.l != f2.l:
			return 0
		ff = GaussFunctionNormed((f1.a + f2.a) / 2, (f1.l + f2.l) // 2)
		return (f1.norm * f2.norm) / ff.norm**2

class GaussFunctionContracted:
	def __init__(self):
		self.fs = []

	def __call__(self, r):
		res = 0
		for c, f in self.fs:
			res += f(r) * c
		return res

	def overlap(self, other):
		res = 0
		for c1, f1 in self.fs:
			for c2, f2 in other.fs:
				fnorm =  GaussFunctionNormed.overlap(f1, f2)
				res += c1 * c2 * fnorm
		return res

	def norm(self):
		return abs(self.overlap(self))**0.5

	def normalize(self):
		n = self.norm()
		self.fs = [(c / n, f) for c, f in self.fs]
		if not [(c,f) for (c,f) in self.fs if c > 0]:
			self.fs = [(-c,f) for (c,f) in self.fs]

	def get_cutoff(self, prec=0.001):
		return self.to_numeric(NumericOperations.loggrid(1e-8, 2e2, 1000)).get_cutoff(prec)

	def to_numeric(self, grid):
		res = NumericFunction()
		for r in grid:
			s = 0
			for c, f in self.fs:
				s += f(r) * c
			res.data.append((r, s))
		res.l = self.fs[0][1].l
		return res

class GaussBasis:
	def __init__(self):
		self.components = {}

	def add_function(self, cg):
		l = cg.fs[0].l
		if l not in list(self.components.keys()):
			self.components[l] = []
		self.components[l].append(cg)
		

class GaussPseudoPotential:
	def __init__(self):
		self.components = {}










	






















