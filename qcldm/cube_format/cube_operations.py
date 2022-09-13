import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..util.elements import ELEMENTS
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
import fortranformat as ff
import numpy as np
from .gaussian_cube import GaussianCube

def rescale_simple(source, target):
	assert not target.is_periodic()
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Rescaling cube')
	logging.info('*********************************************')
	logging.info('')
	res = GaussianCube()
	res.size = [x for x in target.size]
	res.origin = target.origin.copy()
	res.vectors = [x.copy() for x in target.vectors]
	res.celltype = target.celltype
	res.data = np.empty([res.size[x] for x in [0,1,2]])
	res.atoms = target.atoms #TODO: copy
	
	i = 0
	
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				cuboid = target.voxel_cuboid([tx, ty, tz])
				value = source.point_value(cuboid.center)
				res.data[tx,ty,tz] = value
				i += 1
				if i % (res.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, res.data.size))
	return res

def weighted_mean(p1, w1, p2, w2):
	res = []
	for k in range(3):
		res.append((p1[k] * w1 + p2[k] * w2) / (w1 + w2))
	return res

def rescale_medium(source, target):
	assert not target.is_periodic()
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Rescaling cube')
	logging.info('*********************************************')
	logging.info('')
	res = GaussianCube()
	res.size = [x for x in target.size]
	res.origin = target.origin.copy()
	res.vectors = [x.copy() for x in target.vectors]
	res.celltype = target.celltype
	res.data = np.empty([res.size[x] for x in [0,1,2]])
	res.atoms = target.atoms #TODO: copy
	
	i = 0
	
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				cuboid = target.voxel_cuboid([tx, ty, tz])
				value = source.point_value(cuboid.center)
				for v in cuboid.vertices:
					wv = weighted_mean(cuboid.center, 1, v, 9)
					value += source.point_value(wv)
				value /= 9
				res.data[tx,ty,tz] = value
				i += 1
				if i % (res.data.size / 100) == 0:
					logging.debug('  %d of %d' % (i, res.data.size))
	return res

def subtract(cube1, cube2):
	assert cube1.size == cube2.size
	assert cube1.origin == cube2.origin
	assert cube1.vectors == cube2.vectors
	assert cube1.celltype == cube2.celltype
#	assert cube1.atoms == cube2.atoms
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Substracting cubes')
	logging.info('*********************************************')
	logging.info('')
	
	res = GaussianCube()
	res.size = [x for x in cube1.size]
	res.origin = cube1.origin.copy()
	res.vectors = [x.copy() for x in cube1.vectors]
	res.celltype = cube1.celltype
	res.data = np.empty([cube1.size[x] for x in [0,1,2]])
	res.atoms = cube1.atoms #TODO: copy
	
	i = 0
	
	for tx in range(cube1.size[0]):
		for ty in range(cube1.size[1]):
			for tz in range(cube1.size[2]):
				value = cube1.voxel_value([tx, ty, tz]) - cube2.voxel_value([tx, ty, tz])
				res.data[tx,ty,tz] = value
				i += 1
				if i % (res.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, res.data.size))
	return res
				
def masked(target, f):
	assert not target.is_periodic()
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Applying mask')
	logging.info('*********************************************')
	logging.info('')
	res = GaussianCube()
	res.size = [x for x in target.size]
	res.origin = target.origin.copy()
	res.vectors = [x.copy() for x in target.vectors]
	res.celltype = target.celltype
	res.data = np.empty([res.size[x] for x in [0,1,2]])
	res.atoms = target.atoms #TODO: copy
	
	i = 0
	
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				cuboid = target.voxel_cuboid([tx, ty, tz])
				value = target.voxel_value([tx, ty, tz]) if f(cuboid.center) else 0
				res.data[tx,ty,tz] = value
				i += 1
				if i % (res.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, res.data.size))
	return res

def masked_atomsphere(target, atomnum, r):
	center = target.atoms[atomnum - 1].position()
	f = lambda l: ((l[0] - center.x)**2 + (l[1] - center.y)**2 + (l[2] - center.z)**2)**0.5 < r
	return masked(target, f)

def integrate(target):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Integrating')
	logging.info('*********************************************')
	logging.info('')
	res = 0
	i = 0
	dv = target.voxel_volume()
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				res += target.voxel_value([tx, ty, tz]) * dv
				i += 1
				if i % (target.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, target.data.size))
	return res

def integrate_in_sphere(target, atomnum, r):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Integrating')
	logging.info('*********************************************')
	logging.info('')
	res = 0
	i = 0
	dv = target.voxel_volume()
	center = target.atoms[atomnum - 1].position()
	f = lambda l: ((l[0] - center.x)**2 + (l[1] - center.y)**2 + (l[2] - center.z)**2)**0.5 < r
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				cuboid = target.voxel_cuboid([tx, ty, tz])
				if f(cuboid.center):
					res += abs(target.voxel_value([tx, ty, tz]) * dv)
				i += 1
				if i % (target.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, target.data.size))
	return res

def integrate_in_sphere_range(target, atomnum, r, step):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Integrating')
	logging.info('*********************************************')
	logging.info('')
	i = 0
	dv = target.voxel_volume()
	center = target.atoms[atomnum - 1].position()
	funcs = []
	rs = [1.0 * x * step for x in range(0, int(math.ceil(r / step)))]
	for rr in rs:
		f = lambda l,rr=rr: ((l[0] - center.x)**2 + (l[1] - center.y)**2 + (l[2] - center.z)**2)**0.5 < rr
		funcs.append(f)
	res = [0] * len(funcs)
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				cuboid = target.voxel_cuboid([tx, ty, tz])
				for n, f  in enumerate(funcs):
					if f(cuboid.center):
						res[n] += abs(target.voxel_value([tx, ty, tz]) * dv)
				i += 1
				if i % (target.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, target.data.size))
	return rs, res

def integrate_in_sphere_range_opt(target, atomnum, r, step):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Integrating')
	logging.info('*********************************************')
	logging.info('')
	i = 0
	dv = target.voxel_volume()
	center = target.atoms[atomnum - 1].position()
	rs = [1.0 * x * step for x in range(0, int(math.ceil(r / step)))]
	res = [0] * len(rs)
	ls = list(reversed(list(enumerate(rs))))
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				cc = target.voxel_center([tx, ty, tz])
				value = abs(target.voxel_value([tx, ty, tz]) * dv)
				r = ((cc[0] - center.x)**2 + (cc[1] - center.y)**2 + (cc[2] - center.z)**2)**0.5
				for n, r0 in ls:
					if r < r0:
						res[n] += value
					else:
						break
				i += 1
				if i % (target.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, target.data.size))
	return rs, res

def integrate_in_sphere_range_normed(target, atomnum, r, step):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Integrating')
	logging.info('*********************************************')
	logging.info('')
	i = 0
	dv = target.voxel_volume()
	center = target.atoms[atomnum - 1].position()
	rs = [1.0 * x * step for x in range(0, int(math.ceil(r / step)))]
	res = [0] * len(rs)
	res_vol = [0] * len(rs)
	ls = list(enumerate(rs))
	for tx in range(target.size[0]):
		for ty in range(target.size[1]):
			for tz in range(target.size[2]):
				i += 1
				if i % (target.data.size / 10) == 0:
					logging.debug('  %d of %d' % (i, target.data.size))
				
				cc = target.voxel_center([tx, ty, tz])
				value = abs(target.voxel_value([tx, ty, tz]) * dv)
				r = ((cc[0] - center.x)**2 + (cc[1] - center.y)**2 + (cc[2] - center.z)**2)**0.5
				for n, r0 in ls:
					if r < r0:
						res[n] += value
						res_vol[n] += dv
						break

	return rs, [r/v if v != 0 else 0 for r,v in zip(res, res_vol)]

def transformed_copy(target, cell, symop):
	origin = target.origin
	new_origin = cell.apply_symop_to_point(origin, symop)
	
	new_vectors = []
	for v in target.vectors:
		new_v = cell.apply_symop_to_vector(v, symop)
		new_vectors.append(new_v)
	
	new_atoms = []
	for a in target.atoms:
		na = AtomVector(a.name(), cell.apply_symop_to_point(a.position(), symop), a.data())
		new_atoms.append(na)
	
	res = GaussianCube()
	res.atoms = new_atoms
	res.origin = new_origin
	res.vectors = new_vectors
	res.size = target.size
	res.data = target.data
	res.celltype = target.celltype
	return res

	
	
#def trimmed(target):
#	assert not target.is_periodic()
#	logging.info(u'')
#	logging.info(u'*********************************************')
#	logging.info(u'  Trimming cube')
#	logging.info(u'*********************************************')
#	logging.info(u'')
#	
#	
#	res = GaussianCube()
#	res.size = [x for x in target.size]
#	res.origin = target.origin.copy()
#	res.vectors = [x.copy() for x in target.vectors]
#	res.celltype = target.celltype
#	res.data = np.empty([res.size[x] for x in [0,1,2]])
#	res.atoms = target.atoms #TODO: copy
#	
#	i = 0
#	
#	for tx in range(target.size[0]):
#		for ty in range(target.size[1]):
#			for tz in range(target.size[2]):
#				cuboid = target.voxel_cuboid([tx, ty, tz])
#				value = target.voxel_value([tx, ty, tz]) if f(cuboid.center) else 0
#				res.data[tx,ty,tz] = value
#				i += 1
#				if i % (res.data.size / 10) == 0:
#					logging.debug(u'  %d of %d' % (i, res.data.size))
#	return res


#class OverlapDataHolder:
#	def __init__(self):
#		self.v = 0
#		self.o = 0
#		self.n = 0



#def convex_walk(cube, coords, cuboid, holder):
##	print coords
#	if len(coords) == 3:
#		res = cube.check_coord(cuboid, coords, holder)
##		print res
#		return res
#	vectors = []
#	sub_origin = cube.origin.copy()
#	for i in range(3):
#		if i < len(coords):
#			vectors.append(cube.vectors[i])
#			sub_origin += coords[i] * cube.vectors[i]
#		else:
#			vectors.append(cube.vectors[i] * cube.size[i])
#		
#	c = Cuboid(sub_origin._data, vectors)
#	ci = cuboid_intersection(c, cuboid)
#	if not ci:
#		return False
#	if ci.volume() < 1e-14:
#		return False
#	center_coords = cube.find_nearest(ci.mass_center())
#	cnum = len(coords)
#	for i in range(cnum):
#		assert coords[i] == center_coords[i], (str(coords[i])  + " " + str(center_coords[i]))
#	start = center_coords[cnum]
#	found = False
#	for nc in xrange(start, cube.size[cnum]):
#		n_coords = coords + [nc]
#		nf = convex_walk(cube, n_coords, cuboid, holder)
#		found = found or nf
#		if not nf:
#			break
#	for nc in xrange(start - 1, -1, -1):
#		n_coords = coords + [nc]
#		nf = convex_walk(cube, n_coords, cuboid, holder)
#		found = found or nf
#		if not nf:
#			break
#	return found

#	def check_coord(self, nc, coords, holder):
#			oc_origin = np.copy(self.origin._data)
#			for i in range(3):
#					oc_origin += coords[i] * self.vectors[i]._data
#			oc = Cuboid(list(oc_origin), [list(v._data) for v in self.vectors])
#			ci = cuboid_intersection(nc, oc)
#			if ci is None:
#				return False
#			olp = ci.volume()
#			value = self.data[tuple(coords)]
#			holder.v += value * olp
#			holder.o += olp
#			holder.n += 1
#			return olp != 0

#	def rescale(self, n_origin, n_size, n_vectors):
#		logging.info(u'')
#		logging.info(u'*********************************************')
#		logging.info(u'  Rescaling cube')
#		logging.info(u'*********************************************')
#		logging.info(u'')
#		n_data = np.empty([n_size[k] for k in [0,1,2]])
#		i = 0
#		for nx in xrange(n_size[0]):
#			for ny in xrange(n_size[1]):
#				for nz in xrange(n_size[2]):
#					nc_origin = np.copy(n_origin._data)
#					for k in range(3):
#						nc_origin += [nx, ny, nz][k] * n_vectors[k]._data
#					nc = Cuboid(list(nc_origin), [list(v._data) for v in n_vectors])
#					holder = OverlapDataHolder()
#					convex_walk(self, [], nc, holder)
#					if holder.o > 1e-14:
#						n_data[nx,ny,nz] = holder.v / holder.o
#					i += 1
#					if i % (n_data.size / 100) == 0:
#						logging.debug(u'  %d of %d' % (i, n_data.size))
#		self.size = n_size
#		self.vectors = n_vectors
#		self.data = n_data
#		self.origin = n_origin
								

			
			
			
			
			
