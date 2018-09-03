import re, os, sys, math, logging
from math3d import Vector
import numpy as np
from scipy.spatial import Delaunay
import itertools

def tetrahedron_volume(a, b, c, d):
	return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6
	
def hull_volume(pts):
	if not validate_points(pts):
		return 0
	pts = [np.array([p.x, p.y, p.z]) for p in pts]
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

	def center(self):
		return self.origin + reduce(lambda a,b: a+b, self.vectors) / 2

	def cutting_planes(self):
		planes = []
		center = self.center()
		end = self.origin + reduce(lambda a,b: a+b, self.vectors)
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[0], self.vectors[1]), center))
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[1], self.vectors[2]), center))
		planes.append(CuttingPlane(plane_2vectors(self.origin, self.vectors[2], self.vectors[0]), center))
		planes.append(CuttingPlane(plane_2vectors(end, self.vectors[0], self.vectors[1]), center))
		planes.append(CuttingPlane(plane_2vectors(end, self.vectors[1], self.vectors[2]), center))
		planes.append(CuttingPlane(plane_2vectors(end, self.vectors[2], self.vectors[0]), center))
		return planes

def cuboid_intersection(c1, c2):
	pc = PlaneCut()
	for cp in c1.cutting_planes() + c2.cutting_planes():
		pc.add(cp)
	return hull_volume(pc.vertices)


















			
