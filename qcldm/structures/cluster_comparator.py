#!/usr/bin/python
import re, math, logging

from math3d import Vector, Orientation

from atom_vector import AtomVector

def composition(c):
	res = {}
	for a in c:
		if a.name() not in res.keys():
			res[a.name()] = 0
		res[a.name()] += 1
	return res

def lengths(a, c):
	res = []
	for a1 in c:
		res.append(a.distance(a1))
	return sorted(res)

def compare_lenghts(list1, list2, precision):
	for l1, l2 in zip(list1, list2):
		if not rough_equal(l1, l2, precision):
			return False
	return True

def rough_equal(a, b, precision):
	return abs(a - b) < precision

def real_cluster(a, c):
	cr = []
	for a1 in c:
		cr.append(a1.real_atom().relative(a.real_atom()))
	ar = a.real_atom().relative(a.real_atom())
	return ar, cr

def rotate_to_axis(ax1, ax2, c2):
	ax_rot = ax1.position().cross(ax2.position())
	if rough_equal(ax_rot.length, 0, 0.000001):
		if rough_equal(ax1.distance(ax2), 0, 0.000001):
			return c2
		else:
			res = []
			for a in c2:
				res.append(AtomVector(a._name, - a._vector))
			return res
	angle = ax2.position().signed_angle_to(ax1.position(), ax_rot)
	rot = Orientation.new_axis_angle(ax_rot, angle)
	res = []
	for a in c2:
		res.append(AtomVector(a._name, rot * a._vector))
	return res

def rotate_along(ax2, ay1, ay2, c2):
	ax_rot = ax2.position()
	r1 = ay1.position() - ay1.position().projection(ax_rot)
	r2 = ay2.position() - ay2.position().projection(ax_rot)
	angle = r2.signed_angle_to(r1, ax_rot)
	rot = Orientation.new_axis_angle(ax_rot, angle)
	res = []
	for a in c2:
		res.append(AtomVector(a._name, rot * a._vector))
	return res

def reflect(c, ax1, ax2):
	ax_refl = ax1.position().cross(ax2.position())
	res = []
	for a in c:
		res.append(AtomVector(a._name, a._vector - a._vector.projection(ax_refl) * 2))
	return res

def try_compare_two(a1, c1, a2, c2, precision):
	res = True
	res = res and rough_equal(a1.distance(c1[0]), a2.distance(c2[0]), precision)
	res = res and rough_equal(a1.distance(c1[1]), a2.distance(c2[1]), precision)
	res = res and rough_equal(c1[0].distance(c1[1]), c2[0].distance(c2[1]), precision)
	return res

def compare_two(a1, c1, a2, c2, precision):
	return try_compare_two(a1, c1, a2, c2, precision) or try_compare_two(a1, c1, a2, c2[::-1], precision)

def compare_ready(c1, c2, precision):
	for c in c1:
		found = False
		for cc in c2:
			same = c._name == cc._name and rough_equal(c.distance(cc), 0, precision)
			if same:
				found = True
				break
		if not found:
			return False
	return True


def compare_clusters(a1, c1, a2, c2, precision=0.001):
	if a1.name() != a2.name():
		return False

	if len(c1) != len (c2):
		return False

	if len(c1) == 0:
		return True

	if composition(c1) != composition(c2):
		return False

	if not compare_lenghts(lengths(a1, c1), lengths(a2, c2), precision):
		return False

	if len(c1) == 1:
		return True

	if len(c1) == 2:
		return compare_two(a1, c1, a2, c2, precision)

	a1, c1 = real_cluster(a1, c1)
	a2, c2 = real_cluster(a2, c2)

	ax1 = c1[0]
	r = ax1.position().length
	ay1 = c1[1]
	for i in range(1, len(c1)):
		ay1 = c1[i]
		if rough_equal(ax1.position().angle(ay1.position()), math.pi, 0.00001) or rough_equal(ax1.position().angle(ay1.position()), 0, 0.00001):
			continue
		else:
			break
	rr = ay1.position().length
	a = ax1.position().angle(ay1.position())

	for i in range(len(c2)):
		ax2 = c2[i]
		if rough_equal(ax2.position().length, r, precision) and ax1._name == ax2._name:
			c2rot = rotate_to_axis(ax1, ax2, c2)
			ax2 = c2rot[i]
			for j in range(len(c2rot)):
				if j != i:
					ay2 = c2rot[j]
					if rough_equal(ay2.position().length, rr, precision) and rough_equal(ax2.position().angle(ay2.position()), a, precision) and ay1._name == ay2._name:
						c2fin = rotate_along(ax2, ay1, ay2, c2rot)
						if compare_ready(c1, c2fin, precision):
							return True
						ax2 = c2fin[i]
						ay2 = c2fin[j]
						c2fin = reflect(c2fin, ax2, ay2)
						if compare_ready(c1, c2fin, precision):
							return True
	return False


def compare_clusters2(a1, c1, i1, a2, c2, i2, precision=0.001):
	a1, c1 = real_cluster(a1, c1)
	a2, c2 = real_cluster(a2, c2)
	ax1 = c1[i1]
	r = ax1.position().length
	ay1 = None
	rr = 0
	
	for i in range(len(c1)):
		tmp = c1[i]
		if i != i1 and not rough_equal(ax1.position().angle(tmp.position()), math.pi, 0.00001):
			ay1 = tmp
			break
	rr = ay1.position().length
	a = ax1.position().angle(ay1.position())
		
	ax2 = c2[i2]
	if rough_equal(ax2.position().length, ax1.position().length, precision) and ax1._name == ax2._name:
		c2rot = rotate_to_axis(ax1, ax2, c2)
		ax2 = c2rot[i2]
		for j2 in range(len(c2rot)):
			if j2 != i2:
				ay2 = c2rot[j2]
				if rough_equal(ay2.position().length, rr, precision) and rough_equal(ax2.position().angle(ay2.position()), a, precision) and ay1._name == ay2._name:
					c2fin = rotate_along(ax2, ay1, ay2, c2rot)
					if compare_ready(c1, c2fin, precision):
						return True
					ax2 = c2fin[i2]
					ay2 = c2fin[j2]
					c2fin = reflect(c2fin, ax2, ay2)
					if compare_ready(c1, c2fin, precision):
						return True
	return False
	




