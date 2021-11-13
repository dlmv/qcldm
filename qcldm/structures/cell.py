import math, logging, sys
from math3d import Vector
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
np.set_printoptions(suppress=True)
np.core.arrayprint._line_width = 180
np.set_printoptions(edgeitems=10, linewidth=100000)
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

	def validate_rel(self, rel):
		changed = True
		while changed:
			oldrel = rel
			oldrel = [x - 1 if x > 0.5 else x for x in oldrel]
			oldrel = [x + 1 if x <= -0.499999 else x for x in oldrel]
			changed = list(oldrel) != list(rel)
			rel = oldrel
		return rel

	def relative_atom(self, atom, cryst_cell=False):
		vector_basis = np.array([v._data for v in self.vectors])
		if (cryst_cell):
			crysmat = np.array(self.cryst_mat)
			vector_basis = np.transpose(crysmat).dot(vector_basis)
		rel = np.linalg.inv(np.transpose(vector_basis)).dot(atom.position()._data)
		rel = self.validate_rel(rel)
		return AtomVector(atom.name(), Vector(rel), {})

	def true_relative_atom(self, atom, cryst_cell=False):
		vector_basis = np.array([v._data for v in self.vectors])
		if (cryst_cell):
			crysmat = np.array(self.cryst_mat)
			vector_basis = np.transpose(crysmat).dot(vector_basis)
		rel = np.linalg.inv(np.transpose(vector_basis)).dot(atom.position()._data)
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
		return self.cartesian_atom(self.apply_symop(rel, symop), False)

	def atom_in_cell(self, atom):
		ra = self.relative_atom(atom, False)
		ca = self.cartesian_atom(ra)
		ca._data = atom.data()
		return ca

	def validate_point(self, point):
		vector_basis = np.array([v._data for v in self.vectors])
		rel_data = np.linalg.inv(np.transpose(vector_basis)).dot(point._data)
		rel_data = self.validate_rel(rel_data)
		new_data = np.transpose(vector_basis).dot(rel_data)
		return Vector(new_data)

	def apply_symop_to_point(self, point, symop, validate=False):
		vector_basis = np.array([v._data for v in self.vectors])
		n,inv,rot_mat,translator = symop
		rel_data = np.linalg.inv(np.transpose(vector_basis)).dot(point._data)
		new_rel_data = np.array(rot_mat).dot(rel_data)
		new_rel_data += np.array(translator)
		if validate:
			new_rel_data = self.validate_rel(new_rel_data)
		new_data = np.transpose(vector_basis).dot(new_rel_data)
		return Vector(new_data)
		
	def apply_symop_to_vector(self, point, symop):
		vector_basis = np.array([v._data for v in self.vectors])
		n,inv,rot_mat,translator = symop
		rel_data = np.linalg.inv(np.transpose(vector_basis)).dot(point._data)
		new_rel_data = np.array(rot_mat).dot(rel_data)
		new_data = np.transpose(vector_basis).dot(new_rel_data)
		return Vector(new_data)

	def apply_symop(self, atom, symop):
		n,inv,rot_mat,translator = symop
		new_vec = np.array(rot_mat).dot(atom.position()._data)
		new_vec += np.array(translator)
		new_vec = self.validate_rel(new_vec)
		return AtomVector(atom.name(), Vector(new_vec), {})

	def number_in_cell(self, atom):
#		print '===', atom
#		print self.relative_atom(atom)
		for i, a in enumerate(self.atoms):
#			print self.relative_atom(a)
			if self.relative_atom(atom).distance(self.relative_atom(a)) < 1e-4:
				return i
		return -1

	def coord_kernel_matrix(self):
		matrix = []
		for i, k in enumerate(self.assym_n):
			for symop in self.symops:
				n,inv,rot_mat,translator = symop
				a1 = self.apply_symop_cartesian(self.atoms[k], symop)
				l = self.number_in_cell(a1)
				assert l >= 0, "symop resulted in foreign atom"
				if l == k:
					mat = [[0 for xx in range(len(self.assym_n)*3)] for yy in range(3)]
					for row in range(3):
						for col in range(3):
							mat[row][i * 3 + col] += rot_mat[row][col]
						mat[row][i * 3 + row] -= 1
					matrix.extend(mat)
		return linalg.null_space(matrix)
	
	def find_assym(self, atom):
		n = self.number_in_cell(atom)
		assert n >= 0
		for an in self.assym_n:
			for symop in self.symops:
				atom1 = self.apply_symop_cartesian(self.atoms[an], symop)
				l = self.number_in_cell(atom1)
				if l == n:
					return an, symop
		return -1, None

	def coord_basis_matrix(self, atoms):
		matrix = []
		vector_mat = np.transpose([v._data for v in self.vectors])
		for atom in atoms:
			an, symop = self.find_assym(atom)
			n,inv,rot_mat,translator = symop
			assert an >= 0
			mat = [[0 for xx in range(len(self.assym_n)*3)] for yy in range(3)]
			for row in range(3):
				for col in range(3):
					mat[row][self.assym_n.index(an) * 3 + col] += vector_mat.dot(rot_mat)[row][col]
			matrix.extend(mat)
		return np.array(matrix)


	def vector_kernel_matrix(self):
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

	def vector_basis_matrix(self, atoms):
		matrix = []
		for i in range(3 * len(atoms)):
			matrix.append([0] * 9)
		for n, atom in enumerate(atoms):
			ra = self.true_relative_atom(atom)
			for i in range(3):
				for j in range(3):
					matrix[3 * n + i][3 * i + j] += ra.position()._data[j]
		return np.array(matrix)
			
	def cut_by_coords(self, x, atoms):
		k2b = self.coord_kernel_matrix()
		b2x = self.coord_basis_matrix(atoms)
		k = linalg.pinv(k2b).dot(linalg.pinv(b2x).dot(x))
		kb = k2b.dot(k)
		xs = b2x.dot(kb)
		return xs, kb

	def cut_and_expand_by_coords(self, x, atoms, atoms1):
		k2b = self.coord_kernel_matrix()
		b2x = self.coord_basis_matrix(atoms)
		b2x1 = self.coord_basis_matrix(atoms1)
		k = linalg.pinv(k2b).dot(linalg.pinv(b2x).dot(x))
		kb = k2b.dot(k)
		xs1 = b2x1.dot(kb)
		return xs1, kb

	def cut_by_vectors(self, x, atoms):
		k2b = self.vector_kernel_matrix()
		b2x = self.vector_basis_matrix(atoms)
		k = linalg.pinv(k2b).dot(linalg.pinv(b2x).dot(x))
		kb = k2b.dot(k)
		xs = b2x.dot(kb)
		return xs, kb

	def cut_and_expand_by_vectors(self, x, atoms, atoms1):
		k2b = self.vector_kernel_matrix()
		b2x = self.vector_basis_matrix(atoms)
		b2x1 = self.vector_basis_matrix(atoms1)
		k = linalg.pinv(k2b).dot(linalg.pinv(b2x).dot(x))
		kb = k2b.dot(k)
		xs1 = b2x1.dot(kb)
		return xs1, kb

	def modify_by_coords(self, x, atoms):
		x, bx = list(self.cut_by_coords(x, atoms))
		b2x = self.coord_basis_matrix(atoms)
		b2a = self.coord_basis_matrix(self.atoms)
		ax = b2a.dot(bx)
		for i, atom in enumerate(self.atoms):
			v = Vector(np.array(ax[i*3:i*3+3]) + atom.position()._data)
			atom.set_position(v)

	def modify_by_vectors(self, x, atoms):
		x, bx = list(self.cut_by_vectors(x, atoms))
		b2a = self.vector_basis_matrix(self.atoms) # maybe make cell_atoms relative??
		ax = b2a.dot(bx)
		for i, atom in enumerate(self.atoms):
			v = Vector(np.array(ax[i*3:i*3+3]) + atom.position()._data)
			atom.set_position(v)
		for i in range(3):
			for j in range(3):
				self.vectors[j]._data[i] += bx[i * 3 + j]

		
				
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















