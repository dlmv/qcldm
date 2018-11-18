import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..util.elements import ELEMENTS
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
import fortranformat as ff
import numpy as np
from polyhedra_intersection import Cuboid, cuboid_intersection

class OverlapDataHolder:
	def __init__(self):
		self.v = 0
		self.o = 0
		self.n = 0



def convex_walk(cube, coords, cuboid, holder):
#	print coords
	if len(coords) == 3:
		res = cube.check_coord(cuboid, coords, holder)
#		print res
		return res
	vectors = []
	sub_origin = cube.origin.copy()
	for i in range(3):
		if i < len(coords):
			vectors.append(cube.vectors[i])
			sub_origin += coords[i] * cube.vectors[i]
		else:
			vectors.append(cube.vectors[i] * cube.size[i])
		
	c = Cuboid(sub_origin._data, vectors)
	ci = cuboid_intersection(c, cuboid)
	if not ci:
		return False
	if ci.volume() < 1e-14:
		return False
	center_coords = cube.find_nearest(ci.mass_center())
	cnum = len(coords)
	for i in range(cnum):
		assert coords[i] == center_coords[i], (str(coords[i])  + " " + str(center_coords[i]))
	start = center_coords[cnum]
	found = False
	for nc in xrange(start, cube.size[cnum]):
		n_coords = coords + [nc]
		nf = convex_walk(cube, n_coords, cuboid, holder)
		found = found or nf
		if not nf:
			break
	for nc in xrange(start - 1, -1, -1):
		n_coords = coords + [nc]
		nf = convex_walk(cube, n_coords, cuboid, holder)
		found = found or nf
		if not nf:
			break
	return found

class GaussianCube:

	def __init__(self):
		self.atoms = []
		self.origin = Vector(0,0,0)
		self.vectors = []
		self.size = []
		self.data = None

	@staticmethod
	def from_file(name):
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Reading cube from %s' % name)
		logging.info(u'*********************************************')
		logging.info(u'')
		with open(name) as f:
			gc = GaussianCube()
			n = 0
			i = 0
			nat = 9999999999
			lines = f.xreadlines()
			for line in lines:
				if n == 2:
					ls = line.split()
					nat = int(ls[0])
					logging.debug(u'  atoms count: %d' % nat)
					gc.origin = Vector([float(x) for x in ls[1:4]])
				elif 3 <= n < 6:
					ls = line.split()
					nv = int(ls[0])
					vv = Vector([float(x) for x in ls[1:]])
					gc.vectors.append(vv)
					gc.size.append(nv)
				elif 6 <= n < 6+nat:
					ls = line.split()
					name = ELEMENTS[int(ls[0])].symbol
					val = int(float(ls[1]))
					v = Vector([float(x) for x in ls[2:]])
					a = AtomVector(name, v)
					a.data()[AtomKeys.FULL_VALENCE] = val
					a.data()[AtomKeys.ESTIMATED_VALENCE] = Shells.estimate_valence_byname(a.name())
					gc.atoms.append(a)
					if n == 5 + nat:
						gc.data = np.empty([gc.size[x] for x in [0,1,2]])
				elif 6+nat <= n:
					ls = [float(k) for k in line.split()]
					for f in ls:
						x = i / gc.size[1] / gc.size[2]
						y = i / gc.size[1] % gc.size[2]
						z = i % gc.size[1] % gc.size[2]
						gc.data[x,y,z] = f
						i += 1
						if i % (gc.data.size / 10) == 0:
							logging.debug(u'  %d of %d' % (i, gc.data.size))	
				n += 1
			return gc

	def to_file(self, name):
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Writing cube to %s' % name)
		logging.info(u'*********************************************')
		logging.info(u'')
		with open(name, 'w') as f:
			f.write("\n\n")
			f.write("%5d%12.6f%12.6f%12.6f\n" % (len(self.atoms), self.origin.x, self.origin.y, self.origin.z))
			for n, v in zip(self.size, self.vectors):
				f.write("%5d%12.6f%12.6f%12.6f\n" % (n, v.x, v.y, v.z))
			for a in self.atoms:
				f.write("%5d%12.6f%12.6f%12.6f%12.6f\n" % (ELEMENTS[a.name()].number, a.data()[AtomKeys.FULL_VALENCE],  a.position().x, a.position().y, a.position().z))
			i = 0
			lf = ff.FortranRecordWriter('(6e13.5)')
			for x in xrange(self.size[0]):
				for y in xrange(self.size[1]):
					for z in xrange(self.size[2]):
						f.write(lf.write([self.data[x,y,z]]))
						i += 1
						if i % (self.data.size / 10) == 0:
							logging.debug(u'  %d of %d' % (i, self.data.size))
						if z % 6 == 5:
							f.write("\n")
					f.write("\n")
			f.write("\n")

	def find_nearest(self, point):
#		def validate(x, m):
#			x = min(x, m - 1)
#			x = max(x, 0)
#			return x
		a = np.array([k._data for k in self.vectors])
		b = (point - self.origin._data)
		coords = [int(x) for x in np.linalg.solve(a, b)]
#		coords = [validate(int(x), s) for x, s in zip(coords, self.size)]
		return coords

	def check_coord(self, nc, coords, holder):
			oc_origin = np.copy(self.origin._data)
			for i in range(3):
					oc_origin += coords[i] * self.vectors[i]._data
			oc = Cuboid(list(oc_origin), [list(v._data) for v in self.vectors])
			ci = cuboid_intersection(nc, oc)
			if ci is None:
				return False
			olp = ci.volume()
			value = self.data[tuple(coords)]
			holder.v += value * olp
			holder.o += olp
			holder.n += 1
			return olp != 0

	def rescale(self, n_origin, n_size, n_vectors):
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Rescaling cube')
		logging.info(u'*********************************************')
		logging.info(u'')
		n_data = np.empty([n_size[x] for x in [0,1,2]])
		i = 0
		for nx in xrange(n_size[0]):
			for ny in xrange(n_size[1]):
				for nz in xrange(n_size[2]):
					nc_origin = np.copy(n_origin._data)
					for k in range(3):
						nc_origin += [nx, ny, nz][k] * n_vectors[k]._data
					nc = Cuboid(list(nc_origin), [list(v._data) for v in n_vectors])
					holder = OverlapDataHolder()
					convex_walk(self, [], nc, holder)
					if holder.o > 1e-14:
						n_data[nx,ny,nz] = holder.v / holder.o
					i += 1
					if i % (n_data.size / 100) == 0:
						logging.debug(u'  %d of %d' % (i, n_data.size))
		self.size = n_size
		self.vectors = n_vectors
		self.data = n_data
		self.origin = n_origin
								
#	def scale(self, sv):
#		opar = self.params
#		npar = []
#		for s, (n, v) in zip(sv, opar):
#			nn = n / s
#			npar.append([nn, v * (1.0 * n / nn)])

#		odata = self.data
#		ndata = np.empty([npar[x][0] for x in [0,1,2]])
#		
#		for x in xrange(npar[0][0]):
#				for y in xrange(npar[1][0]):
#					for z in xrange(npar[2][0]):
#						val = 0
#						i = 0
#						for dx in xrange(sv[0]):
#							xx = x * sv[0] + dx
#							if xx >= opar[0][0]:
#								continue
#							for dy in xrange(sv[1]):
#								yy = y * sv[1] + dy
#								if yy >= opar[1][0]:
#									continue
#								for dz in xrange(sv[2]):
#									zz = z * sv[2] + dz
#									if zz >= opar[2][0]:
#										continue
#									i += 1
#									val += odata[xx,yy,zz]
#						ndata[x,y,z] = val / i
#						
#		self.params = npar
#		self.data = ndata

#	def multiply(self, mv, exclude_border=True):
#		opar = self.params
#		npar = []
#		for m, (n, v) in zip(mv, opar):
#			if exclude_border:
#				npar.append([n * m - m + 1, v])
#			else:
#				npar.append([n * m, v])

#		odata = self.data
#		ndata = np.empty([npar[x][0] for x in [0,1,2]])
#		
#		for x in range(opar[0][0]):
#				for y in xrange(opar[1][0]):
#					for z in range(opar[2][0]):
#						val = 0
#						i = 0
#						for dx in range(mv[0]):
#							if exclude_border:
#								xx = x + dx * (opar[0][0] - 1)
#							else:
#								xx = x + dx * opar[0][0]
#							for dy in range(mv[1]):
#								if exclude_border:
#									yy = y + dy * (opar[1][0] - 1)
#								else:
#									yy = y + dy * opar[1][0]
#								for dz in range(mv[2]):
#									if exclude_border:
#										zz = z + dz * (opar[2][0] - 1)
#									else:
#										zz = z + dz * opar[2][0]
#									ndata[xx,yy,zz] = odata[x,y,z]
#		
#		oat = self.atoms
#		nat = []
#		
#		for dx in range(mv[0]):
#			for dy in range(mv[1]):
#				for dz in range(mv[2]):
#					for a in oat:
#						v = a.position() + opar[0][0] * opar[0][1] * dx + opar[1][0] * opar[1][1] * dy + opar[2][0] * opar[2][1] * dz
#						aa = AtomVector(a.name(), v)
#						for k in a.data().keys():
#							aa.data()[k] = a.data()[k]
#						nat.append(aa)
#					
#		
#									
#									
#		
#		self.params = npar
#		self.data = ndata
#		self.atoms = nat



						
			
			
			
			
			
			
			
			
			
