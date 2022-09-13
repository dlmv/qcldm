import re, os, sys, math, logging
import numpy as np
from scipy.spatial import Delaunay
from functools import reduce

def tetrahedron_volume(a, b, c, d):
	return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6
	
def hull_volume(pts):
	if not validate_points(pts):
		return 0
	pts = [np.array(p) for p in pts]
	try:
		dt = Delaunay(pts)
	except Exception:
		return 0
	tets = dt.points[dt.simplices]
	vol = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1], 
                                tets[:, 2], tets[:, 3]))
	return vol

def validate_points(pts):
#	print pts
	if len(pts) < 4:
		return False
	i = 0
	i1 = -1
	for ti1 in range(len(pts)):
		if ti1 in [i]:
			continue
		if norm3(sub3(pts[i], pts[ti1])) > 5e-14:
			i1 = ti1
			break
#	print i1

	if i1 == -1:
		return False
		
	v1 = sub3(pts[i1], pts[i])

	i2 = -1
	for ti2 in range(len(pts)):
		if ti2 in [i, i1]:
			continue
		v2 = sub3(pts[ti2], pts[i])
		if norm3(cross3(v1, v2)) > 5e-14:
			i2 = ti2
			break

	if i2 == -1:
		return False
	
#	print i2
			
	p = plane_3points(pts[i], pts[i1], pts[i2])
	
	for ti3 in range(len(pts)):
		if ti3 in [i, i1, i2]:
			continue
		if abs(p.distance(pts[ti3])) > 5e-14:
			return True
				
	
#	for i in range(len(pts) - 3):
#		for i1 in range(i, len(pts) - 2):
#			for i2 in range(i1, len(pts) - 1):
#				for i3 in range(i2, len(pts)):
#					v1 = sub3(pts[i1], pts[i])
#					v2 = sub3(pts[i2], pts[i])
#					v3 = sub3(pts[i3], pts[i])
#					if abs(dot3(v1, cross3(v2, v3))) > 1e-14:
#						return True
	return False

#cross_cache = {}
#dot_cache =  {}

def cross3(n1, n2):
	return (n1[1]*n2[2] - n1[2]*n2[1], n1[2]*n2[0] - n1[0]*n2[2], n1[0]*n2[1] - n1[1]*n2[0])

def dot3(n1, n2):
	return n1[0]*n2[0] + n1[1]*n2[1] + n1[2]*n2[2]
	
#def cross3(n1, n2):
#	if n1 in cross_cache:
#		n1cache = cross_cache[n1]
#		if n2 in n1cache:
#			return n1cache[n2]
#		cross = cross3_internal(n1, n2)
#		n1cache[n2] = cross
#		return cross
#	n1cache = {}
#	cross = cross3_internal(n1, n2)
#	n1cache[n2] = cross
#	cross_cache[n1] = n1cache
#	return cross

#def dot3(n1, n2):
#	print len(dot_cache)
#	if n1 in dot_cache:
#		n1cache = dot_cache[n1]
#		if n2 in n1cache:
#			return n1cache[n2]
#		dot = dot3_internal(n1, n2)
#		n1cache[n2] = dot
#		return dot
#	n1cache = {}
#	dot = dot3_internal(n1, n2)
#	n1cache[n2] = dot
#	dot_cache[n1] = n1cache
#	return dot
	
def add3(n1, n2):
	return (n1[0]+n2[0], n1[1]+n2[1], n1[2]+n2[2])

def add3_3(n1, n2, n3):
	return (n1[0]+n2[0]+n3[0], n1[1]+n2[1]+n3[1], n1[2]+n2[2]+n3[2])

def sub3(n1, n2):
	return (n1[0]-n2[0], n1[1]-n2[1], n1[2]-n2[2])
	
def mul3(n, k):
	return (n[0]*k, n[1]*k, n[2]*k)

def norm3(n):
	return (n[0]**2 + n[1]**2 + n[2]**2)**0.5

class Plane:
	def __init__(self, n, d):
		norm = norm3(n)
		self.n = mul3(n, 1/norm)
		self.d = d / norm

	def distance(self, p):
		return dot3(self.n, p) + self.d

	def __repr__(self):
		return '<Plane: {:.5f}x + {:.5f}y + {:.5f}z + {:.5f} = 0>'.format(self.n[0], self.n[1], self.n[2], self.d)

def plane_3points(p1, p2, p3):
	v1 = sub3(p3, p1)
	v2 = sub3(p2, p1)
	cp = cross3(v1, v2)
	d = -dot3(cp, p3)
	return Plane(cp, d)

def plane_2vectors(o, v1, v2):
	return plane_3points(o, add3(o, v1), add3(o, v2))

class CuttingPlane:
	def __init__(self, plane, inner):
		self.plane = plane
		self.inner = inner

def plane3_intersection(p1, p2, p3):
	div = dot3(p1.n, cross3(p2.n, p3.n))
	if div == 0:
		return None
	return add3_3(mul3(cross3(p2.n, p3.n), -p1.d/div), mul3(cross3(p3.n, p1.n), -p2.d/div), mul3(cross3(p1.n, p2.n), -p3.d/div))

def same_side(p1, p2, plane):
	return plane.distance(p1) * plane.distance(p2) >= 0

def same_side_strict(p1, p2, plane):
	return plane.distance(p1) * plane.distance(p2) > 0

class PlaneCut:
	def __init__(self):
		self.planes = []
		self.vertices = []

	def add(self, plane):
		new_vertices = []
		for v in self.vertices:
			if same_side(v, plane.inner, plane.plane):
				new_vertices.append(v)
		self.vertices = new_vertices
		if len(self.planes) > 1:
			for i in range(len(self.planes) - 1):
				plane1 = self.planes[i]
				for j in range(i + 1, len(self.planes)):
					plane2 = self.planes[j]
					p = plane3_intersection(plane.plane, plane1.plane, plane2.plane)
					if p is None or p in self.vertices:
						continue
					ok = True
					for k in list(range(i)) + list(range(i + 1, j)) + list(range(j + 1, len(self.planes))):
						plane3 = self.planes[k]
						if not same_side(p, plane3.inner, plane3.plane):
							ok = False
							break
					if ok:
						self.vertices.append(p)
		self.planes.append(plane)
		
	def has_vertix(self, point):
		for v in self.vertices:
			same = True
			for i in range(3):
				if abs(v[i] - point[i]) > 1e-10:
					same = False
					break
			if same:
				return True
		return False

	def load_optimized(self, cuboid1, cuboid2):
		self.vertices = []
		self.planes = []
		planes1 = cuboid1.cutting_planes
		planes2 = cuboid2.cutting_planes
		for v in cuboid1.vertices:
			cut = False
			for plane in planes2:
				if not same_side(v, plane.inner, plane.plane):
					cut = True
					break
			if not cut:
				self.vertices.append(v)
		for v in cuboid2.vertices:
			if self.has_vertix(v):
				continue
			cut = False
			for plane in planes1:
				if not same_side(v, plane.inner, plane.plane):
					cut = True
					break
			if not cut:
				self.vertices.append(v)

#		self.planes = planes1 + planes2
				
		for i1 in range(len(planes1)):
			plane1 = planes1[i1]
			for i2 in range(i1, len(planes1)):
				if i2 == i1 + 3:
					continue
				plane2 = planes1[i2]
				for j in range(len(planes2)):
					plane3 = planes2[j]
					p = plane3_intersection(plane1.plane, plane2.plane, plane3.plane)
					if p is None or self.has_vertix(p):
						continue
					ok = True
					for it in range(len(planes1)):
						if it == i1 or it == i2:
							continue
						plane = planes1[it]
						if not same_side(p, plane.inner, plane.plane):
							ok = False
							break
					if ok:
						for jt in range(len(planes2)):
							if jt == j:
								continue
							plane = planes2[jt]
							if not same_side(p, plane.inner, plane.plane):
								ok = False
								break
					if ok:
						self.vertices.append(p)

		for i1 in range(len(planes2)):
			plane1 = planes2[i1]
			for i2 in range(i1, len(planes2)):
				if i2 == i1 + 3:
					continue
				plane2 = planes2[i2]
				for j in range(len(planes1)):
					plane3 = planes1[j]
					p = plane3_intersection(plane1.plane, plane2.plane, plane3.plane)
					if p is None or self.has_vertix(p):
						continue
					ok = True
					for it in range(len(planes2)):
						if it == i1 or it == i2:
							continue
						plane = planes2[it]
						if not same_side(p, plane.inner, plane.plane):
							ok = False
							break
					if ok:
						for jt in range(len(planes1)):
							if jt == j:
								continue
							plane = planes1[jt]
							if not same_side(p, plane.inner, plane.plane):
								ok = False
								break
					if ok:
						self.vertices.append(p)
						
					

	def mass_center(self):
		center = [0,0,0]
		for v in self.vertices:
			for k in range(3):
				center[k] += v[k]
		for k in range(3):
			center[k] /= len(self.vertices)
		return center

	def volume(self):
		return hull_volume(self.vertices)

class Cuboid:
	def __init__(self, origin, vectors):
		self.origin = origin
		self.vectors = vectors
		self.end = add3(self.origin, reduce(lambda a,b: add3(a,b), self.vectors))
		self.center = mul3(add3(self.origin, self.end), 0.5)
		vertices = [tuple(self.origin)]
		for v in self.vectors:
			vertices.append(add3(self.origin, v))
		for v in self.vectors:
			vertices.append(sub3(self.end, v))
		vertices.append(self.end)
		self.vertices = vertices



	@property
	def cutting_planes(self):
		planes = []
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[0], self.vectors[1]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[1], self.vectors[2]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[2], self.vectors[0]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.end, self.vectors[0], self.vectors[1]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.end, self.vectors[1], self.vectors[2]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.end, self.vectors[2], self.vectors[0]), self.center))
		return planes


#def have_possible_intersection(c1, c2):
#	for tc1, tc2 in [[c1, c2], [c2, c1]]:
#		for cp in tc1.cutting_planes():
#			same_found = False
#			for v in tc2.vertices:
#				if same_side_strict(v, cp.inner, cp.plane):
#					same_found = True
#					break
#			if not same_found:
#				return False
#	return True

def have_intersection_chance(c1, c2):
	for i in range(3):
		k1 = [v[i] for v in c1.vertices]
		k2 = [v[i] for v in c2.vertices]
		mink1 = min(k1)
		mink2 = min(k2)
		maxk1 = max(k1)
		maxk2 = max(k2)
		if mink1 > maxk2 or mink2 > maxk1:
			return False
	return True

def cuboid_intersection(c1, c2):
	if not have_intersection_chance(c1, c2):
		return None
#	if not have_possible_intersection(c1, c2):
#		return 0
	pc = PlaneCut()
#	for cp in c1.cutting_planes + c2.cutting_planes:
#		pc.add(cp)
#	vol = pc.volume()
	pc.load_optimized(c1,c2)
#	if abs(pc.volume() - vol) > 1e-8:
#		print pc.volume(), vol
	return pc
















			
