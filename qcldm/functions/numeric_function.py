import re, sys, math
from numpy import logspace, sign

class NumericFunction():
	def __init__(self):
		self.n = 0
		self.l = 0
		self.j = 0.5
		self.data = []

	def grid(self):
		return [r for r, v in self.data]

	def norm(self):
		return (self | self) ** 0.5


	def normalize(self):
		norm = self.norm()
		self.data = [(r, f / norm) for r, f in self.data]

	def num_nodes(self):
		vv = self.data[0][1]
		nn = 0
		for r, v in self.data[1:]:
			if vv * v < 0:
				nn += 1
			if v != 0:
				vv = v
		return nn

	def rescaled(self, grid):
		res = NumericFunction()
		res.n, res.l, res.j = self.n, self.l, self.j
		j = 0
		for i, r in enumerate(grid):
			while j < len(self.data)-1 and self.data[j][0] < r:
				j += 1
			if j == 0:
				res.data.append((r, self.data[0][1]))
			elif r > self.data[j][0]:
				res.data.append((r, 0))
			else:
				r0 = self.data[j-1][0]
				r1 = self.data[j][0]
				vv = (self.data[j-1][1] * (r1-r) + self.data[j][1] * (r-r0)) / (r1 - r0)
				res.data.append((r, vv))
		return res

	def get_cutoff(self, prec=0.001):
		s = 0
		norm = self | self
		for n in reversed(range(len(self.data) - 1)):
			r1 = self.data[n][0]
			r2 = self.data[n+1][0]
			dr = r2 - r1
			r = (r2 + r1) / 2
			fr = ((self.data[n+1][1] * r2) ** 2  + (self.data[n][1] * r1) ** 2) / 2
			c = fr * dr
			s += c
			if abs(s / norm) > prec:
				return r2
		return self.data[-1][0]

	def __abs__(self):
		res = NumericFunction()
		res.n, res.l, res.j = self.n, self.l, self.j
		s = sign(NumericOperations.integrate(self.data))
		res.data = [(r, f * s) for r, f in self.data]
		return res

	def __add__(self, other):
		if isinstance(other, self.__class__):
			assert self.grid() == other.grid(), "Different grids!"
			res = NumericFunction()
			res.n, res.l, res.j = self.n, self.l, self.j
			res.data = [(r, f1 + f2) for (r, f1), (_, f2) in zip(self.data,  other.data)]
			return res
		else:
			raise TypeError("Unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))

	def __sub__(self, other):
		if isinstance(other, self.__class__):
			assert self.grid() == other.grid(), "Different grids!"
			res = NumericFunction()
			res.n, res.l, res.j = self.n, self.l, self.j
			res.data = [(r, f1 - f2) for (r, f1), (_, f2) in zip(self.data,  other.data)]
			return res
		else:
			raise TypeError("Unsupported operand type(s) for -: '{}' and '{}'".format(self.__class__, type(other)))

	def __mul__(self, other):
		if isinstance(other, self.__class__):
			assert self.grid() == other.grid(), "Different grids!"
			res = NumericFunction()
			res.n, res.l, res.j = self.n, self.l, self.j
			res.data = [(r, f1 * f2) for (r, f1), (_, f2) in zip(self.data,  other.data)]
			return res
		elif isinstance(other, (int, float)):
			res = NumericFunction()
			res.n, res.l, res.j = self.n, self.l, self.j
			res.data = [(r, f1 * other) for r, f1 in self.data]
			return res
		else:
			raise TypeError("Unsupported operand type(s) for *: '{}' and '{}'".format(self.__class__, type(other)))

	__rmul__ = __mul__

	def __div__(self, other):
		if isinstance(other, self.__class__):
			raise NotImplementedError("No function division yet")
		if isinstance(other, (int, float)):
			res = NumericFunction()
			res.n, res.l, res.j = self.n, self.l, self.j
			res.data = [(r, f1 / other) for r, f1 in self.data]
			return res
		else:
			raise TypeError("Unsupported operand type(s) for /: '{}' and '{}'".format(self.__class__, type(other)))

	def __neg__(self):
		res = NumericFunction()
		res.n, res.l, res.j = self.n, self.l, self.j
		res.data = [(r, -f1) for r, f1 in self.data]
		return res

	def __or__(self, other):	# <self | other>
		if isinstance(other, self.__class__):
			assert self.grid() == other.grid(), "Different grids!"
			data = [(r, f1 * f2) for (r, f1), (_, f2) in zip(self.data,  other.data)]
			return NumericOperations.integrate(data)
		else:
			raise TypeError("Unsupported operand type(s) for |: '{}' and '{}'".format(self.__class__, type(other)))

	def __rshift__(self, other):	# self >> other = project
		if isinstance(other, self.__class__):
			assert self.grid() == other.grid(), "Different grids!"
			res = (self | other) / (other | other) * other
			return res
		else:
			raise TypeError("Unsupported operand type(s) for >>: '{}' and '{}'".format(self.__class__, type(other)))


class NumericOperations:

	@staticmethod
	def integrate(data):
		res = reduce(lambda s, ((r1, f1), (r2, f2)): s + (f1 * r1**2 + f2 * r2**2) * (r2 - r1) / 2, zip(data[:-1], data[1:]), 0)
		return res

	@staticmethod
	def loggrid(rmin, rmax, n=500):
		return list(logspace(math.log(rmin), math.log(rmax), n, True, math.exp(1)))

	@staticmethod
	def gramm_schmidt(fs):
		res = []
		res.append(fs[0])
		for m in range(1, len(fs)):
			tmp = fs[m]
			for n in range(m):
				tmp = tmp - (tmp >> res[n])
			res.append(tmp)
		for r, f in zip(res, fs):
			r.l, r.n, r.j = f.l, f.n, f.j
		return res

	@staticmethod
	def j_to_ae(p1, p2):
		assert p1.grid() == p2.grid(), "Different grids!"
		if p1.j < p2.j:
			tmp = p1
			p1 = p2
			p2 = tmp
		assert p1.l == p2.l,  'wrong l!'
		assert p1.n == p2.n,  'wrong n!'
		assert p1.j == p1.l + 0.5, 'wrong j1'
		assert p2.j == p2.l - 0.5, 'wrong j1'
		arep = (p2 * p1.l + p1 * (p1.l + 1)) / (2 * p1.l + 1)
		arep.l = p1.l
		arep.n = p1.n
		esop = p1 - p2
		esop.l = p1.l
		esop.n = p1.n
		return arep, esop

	@staticmethod
	def ae_to_j(arep, esop):
		assert arep.grid() == esop.grid(), "Different grids!"
		l = arep.l
		f1 = arep + ((1. * l / (2*l + 1)) * esop)
		f1.l = arep.l
		f1.n = arep.n
		f1.j = arep.l + 0.5
		f2 = arep + ((1. * (-l - 1) / (2*l + 1)) * esop)
		f2.l = arep.l
		f2.n = arep.n
		f2.j = arep.l - 0.5
		return f1, f2

	@staticmethod
	def resort_lnj(ps):
		res = {}
		for p in ps:
			if p.l not in res.keys():
				res[p.l] = {}
			if p.n not in res[p.l].keys():
				res[p.l][p.n] = {}
			if p.j in res[p.l][p.n].keys():
				print 'same values!!!', p.n, p.l, p.j
			res[p.l][p.n][p.j] = p
		return res

	@staticmethod
	def resort_ln(ps):
		res = {}
		for p in ps:
			if p.l not in res.keys():
				res[p.l] = {}
			if p.n not in res[p.l].keys():
				res[p.l][p.n] = {}
			res[p.l][p.n][p.j] = p
		for l in res.keys():
			for n in res[l].keys():
				nj = len(res[l][n].keys())
				if nj == 1:
					res[l][n] = res[l][n][res[l][n].keys()[0]]
				elif nj == 2:
					p = NumericOperations.j_to_ae(res[l][n][res[l][n].keys()[0]], res[l][n][res[l][n].keys()[1]])[0]
					res[l][n] = p
				else:
					print "Wrong j number!"
		return res

	@staticmethod
	def parabolic_smooth(data, limit):
		y0 = limit
		x0 = 0.
		k = None
		pp = []
		for (r1, v1), (r2, v2) in zip(data[:-1], data[1:]):
			if abs(v1) > y0 and abs(v2) <= y0:
				x0 = r2
				y0 = v2
				k = (v2 - v1) / (r2 - r1)
				break
		if k == None:
			return data
		a = k / 2 / x0
		b = y0 - a * x0**2
		newdata = []
		for r, v in data:
			if r < x0:
				v = a * r**2 + b
			newdata.append((r, v))
		return newdata

class CustomFunctions:

	@staticmethod
	def numerize(f, grid):
		res = NumericFunction()
		for r in grid:
			res.data.append((r, f(r)))
		return res

	@staticmethod
	def zero(f, grid):
		return CustomFunctions.numerize(lambda r: 0, grid)

	@staticmethod
	def restricted_sine(grid, rc, nodes):
		f = lambda r: (2. / rc)**0.5 * math.sin(math.pi * (nodes + 1) * r / rc) if rc > r else 0
		return CustomFunctions.numerize(f, grid)







