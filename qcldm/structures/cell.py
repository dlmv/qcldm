import math, logging, sys
from math3d import Vector
import numpy as np
#np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(suppress=True)
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
		self.coord_param_mat = None
		self.vector_param_mat = None
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
			vector_basis = np.transpose(crysmat).dot(vector_basis)
		rel = np.linalg.inv(np.transpose(vector_basis)).dot(atom.position()._data)
		rel = [x - 1 if x > 0.5 else x for x in rel]
		rel = [x + 1 if x <= -0.499999 else x for x in rel]
		return AtomVector(atom.name(), Vector(rel), {})

	def cartesian_atom(self, atom, cryst_cell=False):
		vector_basis = np.array([v._data for v in self.vectors])
		if (cryst_cell):
			crysmat = np.array(self.cryst_mat)
			vector_basis = np.transpose(crysmat).dot(vector_basis)
		unrel = np.transpose(vector_basis).dot(atom.position()._data)
		return AtomVector(atom.name(), Vector(unrel), {})

	def apply_symop_cartesian(self, atom, symop):
		rel = self.relative_atom(atom, False)
		return self.cartesian_atom(self.apply_symop(rel), False)

	def apply_symop(self, atom, symop):
		n,inv,rot_mat,translator = symop
		new_vec = np.array(rot_mat).dot(atom.position()._data)
		new_vec += np.array(translator)
		new_vec = [x - 1 if x > 0.5 else x for x in new_vec]
		new_vec = [x + 1 if x <= -0.499999 else x for x in new_vec]
		return AtomVector(atom.name(), Vector(new_vec), {})

	def number_in_cell(self, atom):
		for i, a in enumerate(self.atoms):
			if self.relative_atom(atom).distance(self.relative_atom(a)) < 1e-5:
				return i
		return -1

	def build_coord_parameters(self):
		matrix = []
		for i in range(len(self.atoms) * len(self.symops) * 3):
			matrix.append([0] * len(self.atoms) * 3)
		for i, a1 in enumerate(self.atoms):
			rela1 = self.relative_atom(a1)
			for k, symop in enumerate(self.symops):
				rela2 = self.apply_symop(rela1, symop)
				n,inv,rot_mat,translator = symop
				j = self.number_in_cell(self.cartesian_atom(rela2))
				assert j >= 0, "symop resulted in foreign atom"
				for row in range(3):
					for col in range(3):
						matrix[i * len(self.symops) * 3 + k * 3 + row][i * 3 + col] += rot_mat[row][col]
				for d in range(3):
					matrix[i * len(self.symops) * 3 + k * 3 + d][j * 3 + d] -= 1
		return linalg.null_space(matrix)

	def build_vector_parameters(self):
		matrix = []
		vector_mat = np.transpose([v._data for v in self.vectors])
		for i in range(9 * len(self.symops)):
			matrix.append([0] * 9)
		for s, symop in enumerate(self.symops):
			n,inv,rot_mat,translator = symop
			rot_mat = np.array(rot_mat)
			cartesian_rot_mat = vector_mat.dot(rot_mat.dot(linalg.inv(vector_mat)))
			for i in range(3):
				for j in range(3):
					for k in range(3):
						matrix[s * 9 + i * 3 + j][3 * i + k] += rot_mat[k][j]
						matrix[s * 9 + i * 3 + j][3 * k + j] -= cartesian_rot_mat[i][k]
		crysmat = np.transpose(self.cryst_mat)
		for i in range(2):
			for j in range(i+1,3):
				m_add = [0] * 9
				for k in range(3):
					m_add[3 * j + k] += crysmat[i][k]
				matrix.append(m_add)
		matrix = np.array(matrix)
		return linalg.null_space(matrix)
#		for i in range(len(kernel[0])):
#			tmp = [[0,0,0],[0,0,0],[0,0,0]]
#			for j in range(len(kernel)):
#				tmp[j % 3][j / 3] = kernel[j][i]#transposed!
#			print np.transpose(self.cryst_mat).dot(tmp)
		
	def cut_by_symmetry(self, vector):
		self.coord_param_mat = self.build_coord_parameters()
		c = np.array(vector)
		c_bas = linalg.pinv(self.coord_param_mat).dot(c)
		c_cut = self.coord_param_mat.dot(c_bas)
		return c_cut

	def cut_by_vectors_symmetry(self, vector):
		self.vector_param_mat = self.build_vector_parameters()
		c = np.array(vector)
		c_bas = linalg.pinv(self.vector_param_mat).dot(c)
		c_cut = self.vector_param_mat.dot(c_bas)
		return c_cut

	def get_vector_basis(self):
		basis = []
		for i in range(3):
			for j in range(3):
				basis.append(self.vectors[j, i])
		return basis

	def get_coord_basis(self):
		full_coords = []
		for atom in self.atoms:
			for k in self.relative_atom(atom).position()._data:
				full_coords.append(k)
		return full_coords

	def set_coord_basis(self, n_coords):
		disps = [n - o for n, o in zip(n_coords, self.get_coord_basis())]
		disps = list(self.cut_vector_by_symmetry(disps))
		coords = [n + o for n, o in zip(disps, self.get_coord_basis())]
		for i, atom in enumerate(self.atoms):
			tmp = AtomVector(atom.name(), Vector(coords[i*3:i*3+3]), atom.data())
			atom.set_position(self.cartesian_atom(tmp).position())
				
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















