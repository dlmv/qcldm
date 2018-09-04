import re, os, sys, math, logging
from math3d import Vector
import numpy as np
from scipy.spatial import Delaunay

def tetrahedron_volume(a, b, c, d):
	return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6
	
def hull_volume(pts):
	if not validate_points(pts):
		return 0
	pts = [p._data for p in pts]
	dt = Delaunay(pts)
	tets = dt.points[dt.simplices]
	vol = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1], 
                                tets[:, 2], tets[:, 3]))
	return vol

def validate_points(pts):
	if len(pts) < 4:
		return False
	for i in range(len(pts) - 3):
		for i1 in range(len(pts) - 2):
			for i2 in range(len(pts) - 1):
				for i3 in range(len(pts)):
					v1 = pts[i1] - pts[i]
					v2 = pts[i2] - pts[i]
					v3 = pts[i3] - pts[i]
					if (v1 * (v2.cross(v3))) != 0:
						return True
	return False
						
					

class Plane:
	def __init__(self, n, d):
		self.n = n / n.length
		self.d = d / n.length

	def distance(self, p):
		return self.n * p + self.d

	def __repr__(self):
		return '<Plane: {:.5f}x + {:.5f}y + {:.5f}z + {:.5f} = 0>'.format(self.n.x, self.n.y, self.n.z, self.d)

def plane_3points(p1, p2, p3):
	v1 = p3 - p1
	v2 = p2 - p1
	cp = v1.cross(v2)
	d = - cp * p3
	return Plane(cp, d)

def plane_2vectors(o, v1, v2):
	return plane_3points(o, o + v1, o + v2)

class CuttingPlane:
	def __init__(self, plane, inner):
		self.plane = plane
		self.inner = inner

def plane3_intersection(p1, p2, p3):
	div = (p1.n * p2.n.cross(p3.n))
	if div == 0:
		return None
	return -(p1.d * (p2.n.cross(p3.n)) + p2.d * (p3.n.cross(p1.n)) + p3.d * (p1.n.cross(p2.n))) / div

def same_side(p1, p2, plane):
	return plane.distance(p1) * plane.distance(p2) >= 0

class PlaneCut:
	def __init__(self):
		self.planes = []
		self.vertices = []

	def add(self, plane):
		new_vertices = []
		for v in self.vertices:
			ok = True
			if not same_side(v, plane.inner, plane.plane):
				ok = False
			if ok:
				new_vertices.append(v)
		self.vertices = new_vertices
		if len(self.planes) > 1:
			for i in range(len(self.planes) - 1):
				plane1 = self.planes[i]
				for j in range(i + 1, len(self.planes)):
					plane2 = self.planes[j]
					p = plane3_intersection(plane.plane, plane1.plane, plane2.plane)
					if not p or p in self.vertices:
						continue
					ok = True
					for k in range(i) + range(i + 1, j) + range(j + 1, len(self.planes)):
						plane3 = self.planes[k]
						if not same_side(p, plane3.inner, plane3.plane):
							ok = False
							break
					if ok:
						self.vertices.append(p)
		self.planes.append(plane)

class Cuboid:
	def __init__(self, origin, vectors):
		self.origin = origin
		self.vectors = vectors
		self.end = self.origin + reduce(lambda a,b: a+b, self.vectors)
		self.center = (self.origin + self.end) / 2
		vertices = [self.origin]
		for v in self.vectors:
			vertices.append(self.origin + v)
		for v in self.vectors:
			vertices.append(self.end - v)
		vertices.append(self.end)
		self.vertices = vertices

	def cutting_planes(self):
		planes = []
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[0], self.vectors[1]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[1], self.vectors[2]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[2], self.vectors[0]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.end, self.vectors[0], self.vectors[1]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.end, self.vectors[1], self.vectors[2]), self.center))
		planes.append(CuttingPlane(plane_2vectors(self.end, self.vectors[2], self.vectors[0]), self.center))
		return planes


def have_possible_intersection(c1, c2):
	for tc1, tc2 in [[c1, c2], [c2, c1]]:
		for cp in tc1.cutting_planes():
			same_found = False
			for v in tc2.vertices:
				if same_side(v, cp.inner, cp.plane):
					same_found = True
					break
			if not same_found:
				return False
	return True

def have_intersection_chance(c1, c2):
	k1 = [v.x for v in c1.vertices]
	k2 = [v.x for v in c2.vertices]
	mink1 = min(k1)
	mink2 = min(k2)
	maxk1 = max(k1)
	maxk2 = max(k2)
	if mink1 > maxk2 or mink2 > maxk1:
		return False

	k1 = [v.y for v in c1.vertices]
	k2 = [v.y for v in c2.vertices]
	mink1 = min(k1)
	mink2 = min(k2)
	maxk1 = max(k1)
	maxk2 = max(k2)
	if mink1 > maxk2 or mink2 > maxk1:
		return False

	k1 = [v.y for v in c1.vertices]
	k2 = [v.y for v in c2.vertices]
	mink1 = min(k1)
	mink2 = min(k2)
	maxk1 = max(k1)
	maxk2 = max(k2)
	if mink1 > maxk2 or mink2 > maxk1:
		return False
	return True

def cuboid_intersection(c1, c2):
	if not have_intersection_chance(c1, c2):
		return 0
#	if not have_possible_intersection(c1, c2):
#		return 0
	pc = PlaneCut()
	for cp in c1.cutting_planes() + c2.cutting_planes():
		pc.add(cp)
	return hull_volume(pc.vertices)


















			
