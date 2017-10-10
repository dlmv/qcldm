import os, logging, numpy
from scipy.linalg import block_diag

from ..structures.atom_vector import AtomKeys
from ..atom.shells import Shells
from ..atom.quantum_numbers import ufu, R2C_Matrix, C2J_Matrix

LMS = 0
LJM = 1

def title(n, oar, t):
	if t == LMS:
		for ll, ns in enumerate(oar):
			if n < ns * 2:
				nn = n / Shells.SHELL_POP[ll] + 1
				n1 = n % Shells.SHELL_POP[ll]
				mm = n1 / 2 - ll
				ss = n1 % 2
				ss =  0.5 if ss == 0 else -0.5
				res = "|%d%s,%d,%.1f>" % (nn, Shells.SHELLS[ll], mm, ss)
				return res, ll, ss
			n -= ns * 2
	elif t == LJM:
		for ll, ns in enumerate(oar):
			if n < ns * 2:
				nn = n / Shells.SHELL_POP[ll] + 1
				n1 = n % Shells.SHELL_POP[ll]
				jj = ll + 0.5 if n1 < (ll + 1) * 2 else ll - 0.5
				if jj == ll - 0.5:
					n1 -= (ll + 1) * 2
				mj = jj - n1
				res = "|%d%s,%.1f,%.1f>" % (nn, Shells.SHELLS[ll], jj, mj)
				return res, ll, jj
			n -= ns * 2
	return 'test'

def matrix_key(a, o1, o2):
	return tuple(sorted([(a.tuple_data(), o1 + 1), (a.tuple_data(), o2 + 1)]))

def build_matrix_s(B, oar):
	res = None
	for ll, ns in enumerate(oar):
		for nn in range(ns / (2 * ll + 1)):
			if res is None:
				res = B.lsmatrix(ll)
			else:
				res = block_diag(res, B.lsmatrix(ll))
	return numpy.matrix(res)

def create_complex_matrix(dms, a, ormat):
	oar = a.data()[AtomKeys.ORBITAL_ARRAY]
	on = a.data()[AtomKeys.ORBITAL_COUNT]
	DM = [[0] * on*2 for x in range(on*2)]
	for i in range(on):
		for j in range(on):
			DM[2*i][2*j] = (dms['a']['a']['re'].get(matrix_key(a, i, j)) or 0) + 1j * (dms['a']['a']['im'].get(matrix_key(a, i, j)) or 0)
			DM[2*i + 1][2*j] = (dms['a']['b']['re'].get(matrix_key(a, i, j)) or 0) + 1j * (dms['a']['b']['im'].get(matrix_key(a, i, j)) or 0)
			DM[2*i][2*j + 1] = (dms['a']['b']['re'].get(matrix_key(a, i, j)) or 0) - 1j * (dms['a']['b']['im'].get(matrix_key(a, i, j)) or 0)
			DM[2*i + 1][2*j + 1] = (dms['b']['b']['re'].get(matrix_key(a, i, j)) or 0) + 1j * (dms['b']['b']['im'].get(matrix_key(a, i, j)) or 0)
	DM = numpy.matrix(DM)
	write_complex_matrix('1_complex_raw.mat', DM, oar, LMS)
	Ubasis = build_matrix_s(ormat, oar)
	DM = ufu(Ubasis, DM)
	write_complex_matrix('2_complex_ordered.mat', DM, oar, LMS)
	Ulm = build_matrix_s(R2C_Matrix(), oar)
	DM = ufu(Ulm, DM)
	write_complex_matrix('3_complex_lms.mat', DM, oar, LMS)
	Ulj = build_matrix_s(C2J_Matrix(), oar)
	DM = ufu(Ulj, DM)
	write_complex_matrix('4_complex_ljm.mat', DM, oar, LJM)
	return DM

def write_complex_matrix(filename, dm, oar, t):
	with open(filename, 'w') as f:
		lm_map = {}
		f.write("\n{:20s}".format("Diagonal elements:"))
		for i in range(len(dm)):
			tt, l, m = title(i, oar, t)
			e = dm[i,i]
			if (l, m) in lm_map.keys():
				lm_map[l, m] += e
			else:
				lm_map[l, m] = e
			f.write("    {:16s}".format(tt))
		f.write("\n" + " " * 20)
		for i in range(len(dm)):
			e = dm[i,i]
			f.write("{:20.6f}".format(e))
		f.write("\n" * 2)
	
		f.write(" " * 20)
		for i in range(len(dm)):
			f.write("    {:16s}".format(title(i, oar, t)[0]))
		f.write("\n")
		for i in range(len(dm)):
			f.write("    {:16s}".format(title(i, oar, t)[0]))
			for j in range(len(dm)):
				e = dm[i,j]
				f.write("{:20.6f}".format(e))
			f.write("\n")

		f.write("\n{:20s}".format("Merged occ:"))
		for (l, m) in sorted(lm_map.keys()):
			f.write("    {:16s}".format("|%d,%.1f>" % (l, m)))
		f.write("\n" + " " * 20)
		for (l, m) in sorted(lm_map.keys()):
			f.write("{:20.6f}".format(lm_map[l, m]))
	











