import math, logging
from math3d import Vector
import numpy as np
import scipy.linalg as linalg
from scipy.linalg import null_space

from cell_neighbours import NeighbourCache
from cluster_comparator import compare_clusters
from atom_vector import AtomVector

def compare_cell_atoms(c1, c2):
	if c1.num != c2.num:
		return -1 if c1.num < c2.num else 1
	for s1, s2 in zip(c1.shifts, c2.shifts):
		if s1 != s2:
			return -1 if s1 < s2 else 1
	return 0

class CellAtom:
	def __init__(self, cell, num, shifts):
		self.num = num
		self.cell = cell
		self.shifts = Vector(shifts)
		self.pos = None

	def position(self):
		if not self.pos:
			tmp = self.cell.atoms[self.num - 1].position().copy()
			for i in range(len(self.cell.vectors)):
				tmp += (self.shifts[i] * self.cell.vectors[i])
			self.pos = tmp
		return self.pos

	def name(self):
		return self.cell.atoms[self.num - 1].name()

	def data(self):
		return self.cell.atoms[self.num - 1].data()

	def shifted(self, shifts):
		return CellAtom(self.cell, self.num, self.shifts + shifts)

	def relative(self, other=None):
		other = other or self
		return self.shifted(-other.shifts)

	def distance(self, other):
		return self.position().dist(other.position())

	def angle(self, a1, a2):
		return (a1.position() - self.position()).angle(a2.position() - self.position()) * 180 / math.pi

	def real_atom(self):
		return AtomVector(self.name(), self.position(), self.data())
		
	def tuple_data(self):
		return (self.num, int(self.shifts[0]), int(self.shifts[1]), int(self.shifts[2]))

	def __hash__(self):
		return hash(self.tuple_data())

	def __str__(self):
		return "cell_atom %s%d(cell %d %d %d)" % (self.name(), self.num, self.shifts[0], self.shifts[1], self.shifts[2])

	def __repr__(self):
		return self.__str__()

	def __lt__(self, other):
		return compare_cell_atoms(self, other) < 0
	def __gt__(self, other):
		return compare_cell_atoms(self, other) > 0
	def __eq__(self, other):
		return compare_cell_atoms(self, other) == 0
	def __le__(self, other):
		return compare_cell_atoms(self, other) <= 0
	def __ge__(self, other):
		return compare_cell_atoms(self, other) >= 0
	def __ne__(self, other):
		return compare_cell_atoms(self, other) != 0

class Cell:

	def __init__(self, atoms, vectors, cryst_mat, symops, assym_n):
		self.vectors = vectors
		self.cryst_mat = cryst_mat
		self.symops = symops
		self.assym_n = assym_n
		self.atoms = atoms
		self.cell = self.shift([0,0,0])
		if self.vectors:
			self.bordered_cell = [a.real_atom() for a in self.shift([0,0,0])]
			start = Vector(0,0,0)
			for i in range(2):
				for j in range(2):
					for k in range(2):
						self.bordered_cell.append(AtomVector('x', start + self.vectors[0] * i + self.vectors[1] * j + self.vectors[2] * k, {}))

		
		if self.vectors:
			self.supercell = self.extended_cell(1)
		else:
			self.supercell = self.cell
		self.neighbours = NeighbourCache(self)
#		self.group_atoms()

	def relative_atom(self, atom, cryst_cell=False):
		vector_basis = np.array([v._data for v in self.vectors])
		if (cryst_cell):
			crysmat = np.array(self.cryst_mat)
			vector_basis = crysmat.dot(vector_basis)
		rel = np.linalg.inv(np.transpose(vector_basis)).dot(atom.position()._data)
		rel = [x - 1 if x > 0.5 else x for x in rel]
		rel = [x + 1 if x < -0.5 else x for x in rel]
		return AtomVector(atom.name(), Vector(rel), {})

	def cartesian_atom(self, atom, cryst_cell=False):
		vector_basis = np.array([v._data for v in self.vectors])
		if (cryst_cell):
			crysmat = np.array(self.cryst_mat)
			vector_basis = crysmat.dot(vector_basis)
		unrel = np.transpose(vector_basis).dot(atom.position()._data)
#		for ca in self.supercell:
#			if abs(ca.position()[0] - unrel[0]) < 1e-5 and abs(ca.position()[1] - unrel[1]) < 1e-5 and abs(ca.position()[2] - unrel[2]) < 1e-5:
#				return ca.relative()
		return AtomVector(atom.name(), Vector(unrel), {})

	def apply_symop(self, atom, symop):
		rel = self.relative_atom(atom, False)
		n,inv,rot_mat,translator = symop
		new_vec = np.array(rot_mat).dot(rel.position()._data)
		new_vec += np.array(translator)
		new_vec = [x - 1 if x > 0.5 else x for x in new_vec]
		new_vec = [x + 1 if x < -0.5 else x for x in new_vec]
		return self.cartesian_atom(AtomVector(atom.name(), Vector(new_vec), {}), False)

	def get_sym_vars(self, n):
		atom = self.atoms[self.assym_n[n]]
		atom = self.cartesian_atom(self.relative_atom(atom, False))
		resmat = []
		for symop in self.symops:
			satom = self.apply_symop(atom, symop)
			print atom, satom
			if satom.distance(atom) < 1e-5:
#				print '========='
				n,inv,rot_mat,translator = symop
				rm = np.array(rot_mat)
				m = rm - np.identity(3)
				for row in m:
					nrow = [x if abs(x) > 1e-5 else 0 for x in row]
					resmat.append(list(nrow))
#				print m
#				print linalg.orth(m)
		resmat = np.array(resmat)
		print resmat
		print linalg.null_space(resmat)
			
	def shift(self, shifts):
		assert self.vectors or shifts == [0, 0, 0], "Trying to shift non-periodic cell!"
		atoms = []
		for n in range(len(self.atoms)):
			at1 = CellAtom(self, n + 1, shifts)
			atoms.append(at1)
		return atoms

	def extended_cell(self, size):
		atoms = []

		c_range = range(-size, size+1)
		for ia in c_range:
			for ib in c_range:
				for ic in c_range:
					tatoms = self.shift([ia, ib, ic])
					for a in tatoms:
						atoms.append(a)
		return atoms

	def check_eq_atoms(self, i, j):
		a1 = self.cell[i]
		a2 = self.cell[j]
		c1 = self.neighbours.first_neighbours(a1)
		c2 = self.neighbours.first_neighbours(a2)
		return compare_clusters(a1, c1, a2, c2, 0.0001)

	def check_atom_in_group(self, i, group):
		if i in group:
			return True
		for j in group:
			if self.check_eq_atoms(i, j):
				return True
		return False

	def group_atoms(self):
		logging.info(u'')
		logging.info(u'************************************************************')
		logging.info(u'          Grouping atoms')
		logging.info(u'************************************************************')
		logging.info(u'')
		grouplist = []
		for i in range(len(self.cell)):
			logging.debug(u'Grouping: %d/%d' % (i+1, len(self.cell)))
			igroup = set([i])
			for group in grouplist:
				if self.check_atom_in_group(i, group):
					igroup = igroup | group
					grouplist.remove(group)
			grouplist.append(igroup)

		print grouplist















