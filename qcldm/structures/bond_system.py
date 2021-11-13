import logging
from ..structures.atom_vector import AtomKeys
from ..matrix.pop_analysis import mulliken_overlap
import numpy as np
np.set_printoptions(precision=3)

class Bond:
	def __init__(self, a1, a2):
		self.delta = a2.shifts - a1.shifts
		self.a1 = a1.relative()
		self.a2 = a2.relative()

	def __str__(self):
		return "%s-%s" % (self.a1, self.a2.shifted(self.delta))

	def __repr__(self):
		return self.__str__()

class BondData:
	def __call__(self, a1, a2):
		return 0.0

class MullikenOverlapBondData(BondData):

	def __init__(self, dm, olp):
		self.dm = dm
		self.olp = olp

	def __call__(self, a1, a2):
		return mulliken_overlap(a1, a2, self.dm, self.olp)

class DumbBondData(BondData):


	def __call__(self, a1, a2):
		return 1.0

class PreloadedBondData(BondData):

	def __init__(self, override_map):
		self.override_map = override_map

	def load_bonds(self, cell):
		bonds = []
		for a in cell.cell:
			nbs = cell.neighbours.first_neighbours(a, self.override_map)
			for nb in nbs:
				bond = Bond(a, nb)
				if a > nb:
					bonds.append(bond)
		return bonds

	def __call__(self, a1, a2):
		delta = a2.shifts - a1.shifts
		for bond, x in self.data:
			if bond.a1 == a1.relative() and bond.a2 == a2.relative() and delta == bond.delta:
				return x
			elif bond.a1 == a2.relative() and bond.a2 == a1.relative() and delta == -bond.delta:
				return -x
		return 0

class LinearSystemChargeTransferBondData(PreloadedBondData):

	def __init__(self, cell, key, override_map):
		PreloadedBondData.__init__(self, override_map)
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Estimating CT from linear system')
		logging.info(u'*********************************************')
		logging.info(u'')


		bonds = self.load_bonds(cell)
		
	
		b = []
		a = []
	

		for atom in cell.cell:
			aline = [0] * len(bonds)
			for i in range(len(bonds)):
				bond = bonds[i]
				if bond.a1 == atom:
					aline[i] = 1
				elif bond.a2 == atom:
					aline[i] = -1
			a.append(aline)
			b.append(atom.data()[key])

		A = np.array(a)
		B = np.array(b)

		# tol copied from matrix_rank
		tol = np.linalg.svd(A, compute_uv=False).max(axis=-1, keepdims=True)[0] * max(A.shape[-2:]) * np.finfo(float).eps

		xs, res, rank, s = np.linalg.lstsq(A, B, tol)
		
		#assert rank < len(b)

		self.data = zip(bonds, xs)








