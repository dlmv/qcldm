import math
import numpy, sympy
from sympy.physics.quantum.cg import CG
from sympy import S

from ..util.mathutils import frange, ufu

class Basis_Matrix:
	def lmatrix(self, l):
		return numpy.identity(2*l + 1)

#lms order: p-1a p-1b p0a p0b p1a p1b
	def lsmatrix(self, l):
		lm = self.lmatrix(l)
		res = [[0] * (4*l + 2) for x in range(4*l + 2)]
		for i in range(2*l + 1):
			for j in range(2*l + 1):
				res[2*i][2*j] = lm[i,j]
				res[2*i + 1][2*j + 1] = lm[i,j]
		return numpy.matrix(res)

class Orbital_Order(Basis_Matrix):

	def order(self, l):
		return range(2 * l + 1)

	def lmatrix(self, l):
		mat = [[0] * (2 * l + 1) for x in range(2 * l + 1)]
		for i, j in enumerate(self.order(l)):
			mat[i][j] = 1
		return numpy.matrix(mat)

def complexToReal(l, mr, mc):
	if mr == 0:
		return 1 if mc == 0 else 0
	else:
		coef = 1j / 2**0.5 if mr < 0 else 1 / 2**0.5
		if mc == -abs(mr):
			return coef
		elif mc == abs(mr):
			return coef * (-1)**mr * math.copysign(1, mr)
		else:
			return 0

def realToComplex(l, mc, mr):
	if mc == 0:
		return 1 if mr == 0 else 0
	else:
		coef = 1 / 2**0.5 if mc < 0 else (-1)**mc / 2**0.5
		if mr == abs(mc):
			return coef
		elif mr == -abs(mc):
			return coef * 1j * math.copysign(1, mc)
		else:
			return 0

def complexToJ(l, m, s, j, mj):
	s = int(round(s * 2))
	j = int(round(j * 2))
	mj = int(round(mj * 2))
	cg = CG(S(l), S(m), S(1)/2, S(s)/2, S(j)/2, S(mj)/2)
	return float(cg.doit().evalf())

class C2R_Matrix(Basis_Matrix):
	def lmatrix(self, l):
		U = []
		for mr in range(-l, l+1):
			line = [0] * (2*l+1)
			for mc in range(-l, l+1):
				line[l + mc] = complexToReal(l, mr, mc)
			U.append(line)
		return numpy.matrix(U)

class R2C_Matrix(Basis_Matrix):
	def lmatrix(self, l):
		U = []
		for mc in range(-l, l+1):
			line = [0] * (2*l+1)
			for mr in range(-l, l+1):
				line[l + mr] = realToComplex(l, mc, mr)
			U.append(line)
		return numpy.matrix(U)

class C2J_Matrix(Basis_Matrix):
	def lsmatrix(self, l):
		U = []
		for j in (l+0.5, l-0.5):
			if j > 0:
				for mj in frange(-j, j+0.5, 1):
					line = [0] * (2*l+1) * 2
					for nm, m in enumerate(range(-l, l+1)):
						for ns, s in enumerate((0.5, -0.5)):
							line[nm * 2 + ns] = complexToJ(l, m, s, j, -mj)
					U.append(line)
		return numpy.matrix(U)

	def lmatrix(self, l):
		raise RuntimeError()

class J2C_Matrix(Basis_Matrix):
	def lsmatrix(self, l):
		U = []
		for m in range(-l, l+1):
			for s in (0.5, -0.5):
				line = [0] * (2*l+1) * 2
				for nj, j in enumerate((l+0.5, l-0.5)):
					if j > 0:
						for nmj, mj in enumerate(frange(-j, j+0.5, 1)):
							line[nj * ((l + 1) * 2) + nmj] = complexToJ(l, m, s, j, -mj)
				U.append(line)
		return numpy.matrix(U)

	def lmatrix(self, l):
		raise RuntimeError()

#print numpy.linalg.inv(C2J_Matrix().lsmatrix(1))


#print ufu(C2J_Matrix().lsmatrix(1), numpy.identity(6))









































