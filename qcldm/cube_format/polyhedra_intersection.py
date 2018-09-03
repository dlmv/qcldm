import re, os, sys, math, logging
from math3d import Vector
import numpy as np
from scipy.spatial import Delaunay

def tetrahedron_volume(a, b, c, d):
	return np.abs(np.einsum('ij,ij->i', a-d, np.cross(b-d, c-d))) / 6
	
def hull_volume(pts):
	if len(pts) < 4:
		return 0
	pts = [np.array([p.x, p.y, p.z]) for p in pts]
	dt = Delaunay(pts)
	tets = dt.points[dt.simplices]
	vol = np.sum(tetrahedron_volume(tets[:, 0], tets[:, 1], 
                                tets[:, 2], tets[:, 3]))
	return vol

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
	print  plane.distance(p1), plane.distance(p2)
	return plane.distance(p1) * plane.distance(p2) >= 0

class PlaneCut:
	def __init__(self):
		self.planes = []
		self.vertices = []

	def add(self, plane):
		new_vertices = []
		for v in self.vertices:
			print v
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











			
