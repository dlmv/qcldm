import math
import numpy, sympy
from sympy.physics.quantum.cg import CG
from sympy import S

def frange(x, y, jump):
	while x < y:
		yield x
		x += jump

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
		if mc == -abs(mr):
			return coef
		elif mc == abs(mr):
			return coef * 1j * math.copysign(1, mr)
		else:
			return 0

def complexToJ(l, m, s, j, mj):
	s = int(round(s * 2))
	j = int(round(j * 2))
	mj = int(round(mj * 2))
	cg = CG(S(l), S(m), S(1)/2, S(s)/2, S(j)/2, S(mj)/2)
	return cg.doit().evalf()

#lms order: p-1a p-1b p0a p0b p1a p1b
def r2c_matrix(l):
	U = []
	for mr in range(-l, l+1):
		for s in range(2):
			line = [0] * (2*l+1) * 2
			for mc in range(-l, l+1):
				line[2 * (l + mc) + s] = realToComplex(l, mc, mr)
			U.append(line)
	return numpy.matrix(numpy.array(U))

def c2j_matrix(l):
	U = []
	for j in (l-0.5, l+0.5):
		for mj in frange(-j, j+0.5, 1):
			line = [0] * (2*l+1) * 2
			for m in range(-l, l+1):
				for sn, s in enumerate((0.5, -0.5)):
					line[2 * (l + m) + sn] = complexToJ(l, m, s, j, mj)
			U.append(line)

	return numpy.matrix(numpy.array(U))

def ufu(u, f):
	return u.getH().dot(f).dot(u)




u = r2c_matrix(1)
#print u
#print '=========='
#print ufu(u, numpy.array([ [1,0,0,0,0,0],[0,0,0,0,0,0],[0,0,1,0,0,0],[0,0,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,0] ]))
u = c2j_matrix(1)
print u
print '==============='
print ufu(u, numpy.array([ [1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1] ]))














































