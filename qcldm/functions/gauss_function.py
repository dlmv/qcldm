import re, sys, math
sys.dont_write_bytecode = True

from numeric_function import NumericFunction

def dfac(n):
	return 1 if n < 2 else reduce(lambda x,y: y*x, range(n,1,-2))

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

class GaussFunctionContracted:
	def __init__(self):
		self.fs = []

	def __call__(self, r):
		res = 0
		for c, f in self.fs:
			res += f(r) * c
		return res

	def norm(self):
		res = 0
		for c1, f1 in self.fs:
			for c2, f2 in self.fs:
				ff = GaussFunctionNormed((f1.a + f2.a) / 2, (f1.l + f2.l) / 2)
				fnorm =  (f1.norm * f2.norm) / ff.norm**2
				res += c1 * c2 * fnorm
		return res**0.5

	def normalize(self):
		n = self.norm()
		self.fs = [(c / n, f) for c, f in self.fs]

	def to_numeric(self, grid):
		res = NumericFunction()
		for r in grid:
			s = 0
			for c, f in self.fs:
				s += f(r) * c
			res.data.append((r, s))
		res.l = self.fs[0][1].l
		return res
















	






















