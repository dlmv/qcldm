import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..util.elements import ELEMENTS
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
from .polyhedra_intersection import Cuboid
import fortranformat as ff
import numpy as np

class GaussianCube:

	TYPE_CLUSTER = 0
	TYPE_PERIODIC_FULL = 1
	TYPE_PERIODIC_WITH_BORDER = 2

	def __init__(self):
		self.atoms = []
		self.origin = Vector(0,0,0)
		self.vectors = []
		self.size = []
		self.data = None
		self.celltype = self.TYPE_CLUSTER

	def is_periodic(self):
		return self.celltype != self.TYPE_CLUSTER

	@staticmethod
	def from_file(name):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Reading cube from %s' % name)
		logging.info('*********************************************')
		logging.info('')
		with open(name) as f:
			gc = GaussianCube()
			n = 0
			i = 0
			nat = 9999999999
			lines = f
			for line in lines:
				if n == 2:
					ls = line.split()
					nat = int(ls[0])
					logging.debug('  atoms count: %d' % nat)
					gc.origin = Vector([float(x) for x in ls[1:4]])
				elif 3 <= n < 6:
					ls = line.split()
					nv = int(ls[0])
					vv = Vector([float(x) for x in ls[1:]])
					gc.vectors.append(vv)
					gc.size.append(nv)
				elif 6 <= n < 6+nat:
					ls = line.split()
					num = int(ls[0])
					name = ELEMENTS[num].symbol if num != 0 else 'Q'
					val = float(ls[1])
					v = Vector([float(x) for x in ls[2:]])
					a = AtomVector(name, v)
					a.data()[AtomKeys.FULL_VALENCE] = val
					gc.atoms.append(a)
					if n == 5 + nat:
						gc.data = np.empty([gc.size[x] for x in [0,1,2]])
				elif 6+nat <= n:
					ls = [float(k) for k in line.split()]
					for f in ls:
						x = i / gc.size[2] / gc.size[1]
						y = i / gc.size[2] % gc.size[1]
						z = i % gc.size[2] % gc.size[1]
						gc.data[x,y,z] = f
						i += 1
						if i % (gc.data.size / 10) == 0:
							logging.debug('  %d of %d' % (i, gc.data.size))	
				n += 1
			return gc

		
	def to_file(self, name):
		ranges = []
		for k in range(3):
			ranges.append([0, self.size[k]])
		self.to_file_custom(name, ranges)

	def to_file_custom(self, name, ranges):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Writing cube to %s' % name)
		logging.info('*********************************************')
		logging.info('')
		with open(name, 'w') as f:
			f.write("\n\n")
			f.write("%5d%12.6f%12.6f%12.6f\n" % (len(self.atoms), self.origin.x, self.origin.y, self.origin.z))
			for n, v in zip(self.size, self.vectors):
				f.write("%5d%12.6f%12.6f%12.6f\n" % (n, v.x, v.y, v.z))
			for a in self.atoms:
				f.write("%5d%12.6f%12.6f%12.6f%12.6f\n" % (ELEMENTS[a.name()].number if a.name() != 'Q' else 0, a.data()[AtomKeys.FULL_VALENCE],  a.position().x, a.position().y, a.position().z))
			i = 0
			lf = ff.FortranRecordWriter('(6e13.5)')
			size = (ranges[0][1] - ranges[0][0]) * (ranges[1][1] - ranges[1][0]) * (ranges[2][1] - ranges[2][0])
			for x in range(ranges[0][0], ranges[0][1]):
				for y in range(ranges[1][0], ranges[1][1]):
					for z in range(ranges[2][0], ranges[2][1]):
						f.write(lf.write([self.voxel_value([x,y,z])]))
						i += 1
						if i % (size / 10) == 0:
							logging.debug('  %d of %d' % (i, size))
						if z % 6 == 5:
							f.write("\n")
					f.write("\n")
			f.write("\n")

	def voxel_value(self, coords):
		assert len(coords) == 3
		n_coords = []
		for k in range(3):
			if coords[k] >= 0 and coords[k] < self.size[k]:
				n_coords.append(coords[k])
			elif self.celltype == self.TYPE_CLUSTER:
				return 0
			elif self.celltype == self.TYPE_PERIODIC_FULL:
				n_coords.append(coords[k] % self.size[k])
			elif self.celltype == self.TYPE_PERIODIC_WITH_BORDER:
				n_coords.append(coords[k] % (self.size[k] - 1))
		return self.data[tuple(n_coords)]

	def voxel_cuboid(self, coords):
		assert len(coords) == 3
		n_origin = np.copy(self.origin._data)
		for k in range(3):
				n_origin += coords[k] * self.vectors[k]._data
		return Cuboid(list(n_origin), [list(v._data) for v in self.vectors])
		
	def voxel_center(self, coords):
		assert len(coords) == 3
		res = np.copy(self.origin._data)
		for k in range(3):
				res += (coords[k] + 0.5) * self.vectors[k]._data
		return res
		

	def find_voxel(self, point):
		a = np.array([k._data for k in self.vectors]).transpose()
		b = (point - self.origin._data)
		coords = [int(math.floor(x)) for x in np.linalg.solve(a, b)]
		return coords

	def point_value(self, point):
		return self.voxel_value(self.find_voxel(np.array(point)))

	def voxel_volume(self):
		return self.vectors[0] * (self.vectors[1].cross(self.vectors[2]))
	



						
			
			
			
			
			
			
			
			
			
