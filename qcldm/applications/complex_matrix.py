import os, logging, numpy, scipy
from scipy.linalg import block_diag

from ..structures.atom_vector import AtomKeys
from ..atom.shells import Shells
from ..atom.harmonics import C2R_Matrix, J2C_Matrix
from ..util.mathutils import ufu


LMS = 0
LJM = 1
LMS_OLP = 2

def title(n, oar, t):
	if t == LMS or t == LMS_OLP:
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

def convert_atom_matrix(DM, OLP, a, ormat):
	oar = a.data()[AtomKeys.ORBITAL_ARRAY]
	DM = numpy.matrix(DM)
	OLP = numpy.matrix(OLP)
	write_complex_matrix('1_olp_input.mat', OLP, oar, LMS_OLP)
	write_complex_matrix('1_dm_input.mat', DM, oar, LMS)
	S12 = scipy.linalg.sqrtm(OLP)
	DM = S12.dot(DM).dot(S12)
	write_complex_matrix('2_dm_overlapped.mat', DM, oar, LMS)
	Ubasis = build_matrix_s(ormat, oar)
	DM = ufu(Ubasis, DM)
	write_complex_matrix('3_dm_ordered.mat', DM, oar, LMS)
	Ulm = build_matrix_s(C2R_Matrix(), oar)
	DM = ufu(Ulm, DM)
	write_complex_matrix('4_dm_lms.mat', DM, oar, LMS)
	Ulj = build_matrix_s(J2C_Matrix(), oar)
	DM = ufu(Ulj, DM)
	write_complex_matrix('5_dm_ljm.mat', DM, oar, LJM)
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
		if t == LMS_OLP:
			return
		f.write("\n\nMerged occ:\n")
		s = 0
		for (l, m) in sorted(lm_map.keys()):
			f.write("{:10s}".format("|%s,%.1f>" % (Shells.SHELLS[l], m)))
			f.write("{:20.6f}".format(lm_map[l, m]))
			f.write("\n")
			s += lm_map[l, m]
		f.write("\n\nTotal occ:")
		f.write("{:20.6f}".format(s))
	











