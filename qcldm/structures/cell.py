import math, logging
from math3d import Vector

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
		return "atom %s%d(cell %d %d %d)" % (self.name(), self.num, self.shifts[0], self.shifts[1], self.shifts[2])

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

	def __init__(self, atoms, vectors):
		self.vectors = vectors
		self.atoms = atoms
		self.cell = self.shift([0,0,0])
		if self.vectors:
			self.supercell = self.extended_cell()
		else:
			self.supercell = self.cell
		self.neighbours = NeighbourCache(self)
#		self.group_atoms()

	def shift(self, shifts):
		assert self.vectors or shifts == [0, 0, 0], "Trying to shift non-periodic cell!"
		atoms = []
		for n in range(len(self.atoms)):
			at1 = CellAtom(self, n + 1, shifts)
			atoms.append(at1)
		return atoms

	def extended_cell(self):
		atoms = []
		c_range = (-1, 0, 1)
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















