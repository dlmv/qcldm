#!/usr/bin/python
import math, logging, re

from ..util.elements import ELEMENTS
from ..util.units import Units


def check_distance(a1, a2, r):
	if abs(a1.position().x - a2.position().x) > r:
		return False
	if abs(a1.position().y - a2.position().y) > r:
		return False
	if abs(a1.position().z - a2.position().z) > r:
		return False
	return a1.distance(a2) < r

NONMETALS = ['H', 'C', 'O', 'F', 'P', 'S', 'Cl', 'Se', 'Br', 'I']

class NeighbourCache:
	def __init__(self, c):
		self.cell = c
		self.cache = {}

	def neighbours_cluster(self, centers, layers, override_map):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Cutting cluster for: '  + str(centers))
		logging.info('*********************************************')
		logging.info('')
		shells = [centers]
		for i in range(layers):
			self.expand_neighbours(shells, override_map)
		return shells

	def load_all(self, override_map):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Loading neighbours for the whole cell')
		logging.info('*********************************************')
		logging.info('')
		for i in range(len(self.cell.cell)):
			c = self.cell.cell[i]
			logging.debug('Finding neighbours %d/%d for: %s' % (i+1, len(self.cell.cell), str(c)))
			logging.debug('  Neighbours count: %d' % (len(self.first_neighbours(c, override_map))))

	def write_neighbours(self, name):
		self.load_all()
		with open(name, 'w') as f:
			for c in self.cell.cell:
				nbs = self.first_neighbours(c)
				for nb in nbs:
					f.write('%3d : %3d %3d %3d %3d\n' % (c.num, nb.num, nb.shifts[0], nb.shifts[1], nb.shifts[2]))


	def read_neighbours(self, name):
		self.cache = {}
		with open(name, 'r') as f:
			for line in f.readlines():
				ls  = [_f for _f in re.split(':', line) if _f]
				ls1 = [_f for _f in re.split('\s*', ls[1]) if _f]
				num = int(ls[0])
				if num not in list(self.cache.keys()):
					self.cache[num] = []
				nb = self.cell.cell[int(ls1[0])].shifted(int(ls1[1]), int(ls1[2]), int(ls1[3]))
				self.cache[num].append(nb)

	def expand_neighbours(self, shells, override_map):
		border = shells[-1]
		newborder = []
		for b in border:
			nbs = self.first_neighbours(b, override_map)
			for nb in nbs:
				found = False
				for shell in shells:
					if nb in shell:
						found = True
				if nb in newborder:
					found = True
				if not found:
					newborder.append(nb)
		shells.append(newborder)

	def expand_covalent_border(self, shells, override_map, limit=5):
		i = limit
		while i > 0:
			border = shells[-1]
			newborder = []
			for b in border:
				if b.name() not in NONMETALS:
					continue
				nbs = self.first_neighbours(b, override_map)
				for nb in nbs:
					if nb.name() not in NONMETALS:
						continue
					found = False
					for shell in shells:
						if nb in shell:
							found = True
					if nb in newborder:
						found = True
					if not found:
						newborder.append(nb)
			if not newborder:
				return
			border.extend(newborder)
			i -= 1
		assert False, 'expand limit exceeded'
			

	def merge_from(self, other):
		for n in list(other.cache.keys()):
			if n not in list(self.cache.keys()):
				self.cache[n] = []
			for other_atom in other.cache[n]:
				self_atom = self.cell.cell[other_atom.num - 1].shifted(other_atom.shifts)
				if self_atom not in self.cache[n]:
					self.cache[n].append(self_atom)

	def first_neighbours(self, center, override_map):
		atoms = []
		shifts = center.shifts
		center = center.shifted(-shifts)
		if center.num in list(self.cache.keys()):
			atoms = self.cache[center.num]
		else:
			for a0 in center.cell.supercell:
				if a0 != center:
					if self.test_neighbour(center, a0, override_map):
						atoms.append(a0)
			self.cache[center.num] = atoms
		satoms = []
		for a in atoms:
			satoms.append(a.shifted(shifts))
		return satoms


	def test_neighbour(self, center, a, override_map):
		rc = (ELEMENTS[center.name()].covrad + ELEMENTS[a.name()].covrad) / Units.UNIT
		rc = rc * 0.8 + 0.8
		bond = tuple(sorted([center.name(), a.name()]))
		if bond in list(override_map.keys()):
			rc = override_map[bond]
		if not check_distance(center, a, rc):
			return False
		r = center.distance(a)
		for b in center.cell.supercell:
			if b != center and b != a and check_distance(center, b, r) and check_distance(a, b, r) and b.angle(center, a) > 90:
				return False
		return True



