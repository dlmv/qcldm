import re, sys, math
from numeric_function import NumericFunction
sys.dont_write_bytecode = True

class Potential():
	def __call__(self, f):
		raise NotImplementedError('!!!!')

class SimplePotential(Potential):

	def __init__(self, f):
		self.core = f

class LocalPotential(SimplePotential):

	def __call__(self, f):
		if isinstance(f, NumericFunction):
			return self.core * f
		if isinstance(f, LocalPotential):
			return LocalPotential(self.core * f.core)
		elif isinstance(f, SequencePotential):
			return SequencePotential(f.pps + [self])
		elif isinstance(f, Potential):
			return SequencePotential([f, self])
		else:
			raise TypeError("Unsupported operand type(s) for (): '{}' and '{}'".format(self.__class__, type(other)))

	def __add__(self, other):
		if isinstance(other, LocalPotential):
			return LocalPotential(self.core + other.core)
		elif isinstance(other, SumPotential):
			return SumPotential(other.pps + [self])
		elif isinstance(other, Potential):
			return SumPotential([other, self])
		else:
			raise TypeError("Unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))

	__radd__ = __add__

class ProjectorPotential(SimplePotential):

	def __call__(self, f):
		if isinstance(f, NumericFunction):
			return f >> self.core
		elif isinstance(f, SequencePotential):
			return SequencePotential(f.pps + [self])
		elif isinstance(f, Potential):
			return SequencePotential([f, self])
		else:
			raise TypeError("Unsupported operand type(s) for (): '{}' and '{}'".format(self.__class__, type(other)))

	def __add__(self, other):
		if isinstance(other, SumPotential):
			return SumPotential(other.pps + [self])
		elif isinstance(other, Potential):
			return SumPotential([other, self])
		else:
			raise TypeError("Unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))

	__radd__ = __add__

class ComplexPotential(Potential):
	def __init__(self, pps):
		self.pps = pps

class SequencePotential(ComplexPotential):
	def __call__(self, f):
		if isinstance(f, NumericFunction):
			return reduce(lambda cf, p: p(cf), self.pps, f)
		elif isinstance(f, SequencePotential):
			return SequencePotential(f.pps + self.pps)
		elif isinstance(f, Potential):
			return SequencePotential([f] + self.pps)
		else:
			raise TypeError("Unsupported operand type(s) for (): '{}' and '{}'".format(self.__class__, type(other)))



	def __add__(self, other):
		if isinstance(other, SumPotential):
			return SumPotential(other.pps + [self])
		elif isinstance(other, Potential):
			return SumPotential([other, self])
		else:
			raise TypeError("Unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))

	__radd__ = __add__

class SumPotential(ComplexPotential):
	def __call__(self, f):
		if isinstance(f, NumericFunction):
			return reduce(lambda cf, p: cf + p(f), self.pps, f*0)
		elif isinstance(f, SequencePotential):
			return SequencePotential(f.pps + [self])
		elif isinstance(f, Potential):
			return SequencePotential([f, self])
		else:
			raise TypeError("Unsupported operand type(s) for (): '{}' and '{}'".format(self.__class__, type(other)))

	def __add__(self, other):
		if isinstance(other, SumPotential):
			return SumPotential(other.pps + self.pps)
		elif isinstance(other, Potential):
			return SumPotential([other] + self.pps)
		else:
			raise TypeError("Unsupported operand type(s) for +: '{}' and '{}'".format(self.__class__, type(other)))

	__radd__ = __add__

def grecp_lj_sep(pps, vp):
	res = SumPotential([])

	for cp, cf in pps:
		p = ProjectorPotential(cf)
		u = LocalPotential(cp - vp)
		res += p(u)
		res += u(p)

	for cp1, cf1 in pps:
		for cp2, cf2 in pps:
			p1 = ProjectorPotential(cf1)
			p2 = ProjectorPotential(cf2)
			u = LocalPotential(vp - (0.5 * (cp1 + cp2)))
			res += p1(u(p2))

	return res

def create_pseudo_potential(val, pps, loc, grecp):
	lp = LocalPotential(val - loc)
	return lp + grecp_lj_sep(pps, val) if grecp else lp

def project_v(f, core, v):
	return (core | v(f)) / (core | v(core)) * core

def gramm_schmidt_v(fs, v):
	res = []
	res.append(fs[0])
	for m in range(1, len(fs)):
		tmp = fs[m]
		for n in range(m):
			tmp = tmp - project_v(tmp, res[n], v)
		res.append(tmp)
	return res

def blochl_form(l, fs, p, normalize=False):
	nfs = []
	for pf in fs:
		pf.normalize()
	fs = gramm_schmidt_v(fs, p)
	projectors = []
	proj_energies = []
	for pf in fs:
		pf.normalize()
		e = pf | p(pf)
		norm = (abs(e) ** 0.5) if normalize else 1
		e = math.copysign(1, e)  if normalize else e
		proj_energies.append(1/e)
		proj = p(pf) / norm
		projectors.append(proj)
	return projectors, proj_energies













