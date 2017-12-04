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

class NeighbourCache:
	def __init__(self, c):
		self.cell = c
		self.cache = {}

	def neighbours_cluster(self, num, layers):
		center = self.cell.cell[num - 1]
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Cutting cluster for: '  + str(center))
		logging.info(u'*********************************************')
		logging.info(u'')
		border = [center]
		shells = [border]
		self.load_neighbours_internal(border, layers, shells)
		return shells

	def load_all(self):
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Loading neighbours for the whole cell')
		logging.info(u'*********************************************')
		logging.info(u'')
		for i in range(len(self.cell.cell)):
			c = self.cell.cell[i]
			logging.debug(u'Finding neighbours %d/%d for: %s' % (i+1, len(self.cell.cell), str(c)))
			logging.debug(u'  Neighbours count: %d' % (len(self.first_neighbours(c))))

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
				ls  = filter(None, re.split(':', line))
				ls1 = filter(None, re.split('\s*', ls[1]))
				num = int(ls[0])
				if num not in self.cache.keys():
					self.cache[num] = []
				nb = self.cell.cell[int(ls1[0])].shifted(int(ls1[1]), int(ls1[2]), int(ls1[3]))
				self.cache[num].append(nb)

	def load_neighbours_internal(self, border, layers, shells):
		if layers == 0:
			logging.debug(u'Cluster done!')
			return []
		newborder = []
		logging.debug(u'Finding next cluster shell')
		for b in border:
			nbs = self.first_neighbours(b)
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
		self.load_neighbours_internal(newborder, layers - 1, shells)

	def merge_from(self, other):
		for n in other.cache.keys():
			if n not in self.cache.keys():
				self.cache[n] = []
			for other_atom in other.cache[n]:
				self_atom = self.cell.cell[other_atom.num - 1].shifted(other_atom.shifts)
				if self_atom not in self.cache[n]:
					self.cache[n].append(self_atom)

	def first_neighbours(self, center):
		atoms = []
		shifts = center.shifts
		center = center.shifted(-shifts)
		if center.num in self.cache.keys():
			atoms = self.cache[center.num]
		else:
			for a0 in center.cell.supercell:
				if a0 != center:
					if self.test_neighbour(center, a0):
						atoms.append(a0)
			self.cache[center.num] = atoms
		satoms = []
		for a in atoms:
			satoms.append(a.shifted(shifts))
		return satoms


	def test_neighbour(self, center, a):
		rc = (ELEMENTS[center.name()].covrad + ELEMENTS[a.name()].covrad) / Units.UNIT
		if not check_distance(center, a, rc * 0.8 + 0.8):
			return False
		r = center.distance(a)
		for b in center.cell.supercell:
			if b != center and b != a and check_distance(center, b, r) and check_distance(a, b, r) and b.angle(center, a) > 90:
				return False
		return True



