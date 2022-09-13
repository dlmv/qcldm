#!/usr/bin/python

import sys, os, logging

#from ..openmx_format.dat_format import DAT_INPUT

def corresponding_atom(a, cell2):
	return cell2.cell[a.num - 1].shifted(a.shifts)

def compare_cells(c1, c2):
	if len(c1.cell) != len(c2.cell):
		print('different atom number!!')
		return
	for a1, a2 in zip(c1.cell, c2.cell):
		if a1.name() != a2.name():
			print('different atom types!')
			return
	c1.neighbours.load_all({})
	c2.neighbours.load_all({})
	c1.neighbours.merge_from(c2.neighbours)

	sd = 0
	mxd = 0

	n = 0

	a, nb = 0,0

	for a1 in c1.cell:
		a2 = corresponding_atom(a1, c2)
		for n1 in c1.neighbours.first_neighbours(a1, {}):
			n2 = corresponding_atom(n1, c2)
			r1 = a1.distance(n1)
			r2 = a2.distance(n2)
#			d = ((r1 - r2) / r1)**2
			d = (r1 - r2)**2
			if abs(r1-r2) > 2:
				continue
			if abs(r1-r2) > 0.01:
				print("%7s%6.3f%6.3f" % ("%s%d" % (n1.name(), n1.num), r1, r2))
			sd += d
			if mxd < d and n1 in c1.cell:
				mxd = d
				a, nb = a1, n1
			n += 1
		for n1 in c1.neighbours.first_neighbours(a1, {}):
			n2 = corresponding_atom(n1, c2)
			for m1 in c1.neighbours.first_neighbours(a1, {}):
				if n1 >= m1:
					continue
				m2 = corresponding_atom(m1, c2)
				an1 = a1.angle(n1, m1)
				an2 = a2.angle(n2, m2)
				if abs(an1-an2) > 15:
					continue
				if abs(an1-an2) > 3:
					print("%7s%5s%8.3f%8.3f" % ("%s%d" % (n1.name(), n1.num), "%s%d" % (m1.name(), m1.num), an1, an2))
	sd = (sd / n) ** 0.5
	mxd = mxd ** 0.5
	return sd, mxd, a, nb


def compare_cells3(c1, c2, c3):
	if len(c1.cell) != len(c2.cell) or len(c1.cell) != len(c3.cell):
		print('different atom number!!')
		return
	for a1, a2, a3 in zip(c1.cell, c2.cell, c3.cell):
		if a1.name() != a2.name() or  a1.name() != a3.name():
			print('different atom types!')
			return
	c1.neighbours.load_all({})
	c2.neighbours.load_all({})
	c3.neighbours.load_all({})
	c1.neighbours.merge_from(c2.neighbours)
	c1.neighbours.merge_from(c3.neighbours)

	for a1 in c1.cell:
		a2 = corresponding_atom(a1, c2)
		a3 = corresponding_atom(a1, c3)
		for n1 in c1.neighbours.first_neighbours(a1, {}):
			n2 = corresponding_atom(n1, c2)
			n3 = corresponding_atom(n1, c3)
			r1 = a1.distance(n1)
			r2 = a2.distance(n2)
			r3 = a3.distance(n3)
			if abs(r1-r2) > 0.05 or abs(r1-r3) > 0.05 or abs(r3-r2) > 0.05:
				print("%7s%6.3f%6.3f%6.3f" % ("%s%d" % (n1.name(), n1.num), r1, r2, r3))
		for n1 in c1.neighbours.first_neighbours(a1, {}):
			n2 = corresponding_atom(n1, c2)
			n3 = corresponding_atom(n1, c3)
			for m1 in c1.neighbours.first_neighbours(a1, {}):
				if n1 >= m1:
					continue
				m2 = corresponding_atom(m1, c2)
				m3 = corresponding_atom(m1, c3)
				an1 = a1.angle(n1, m1)
				an2 = a2.angle(n2, m2)
				an3 = a2.angle(n3, m3)
				if abs(an1-an2) > 3 or abs(an1-an3) > 3 or abs(an3-an2) > 3:
					print("%7s%5s%8.3f%8.3f%8.3f" % ("%s%d" % (n1.name(), n1.num), "%s%d" % (m1.name(), m1.num), an1, an2, an3))



def get_dispacements(c1, c2):
	if len(c1.cell) != len(c2.cell):
		print('different atom number!!')
		return
	for a1, a2 in zip(c1.cell, c2.cell):
		if a1.name() != a2.name():
			print('different atom types!')
			return

	sd = 0
	mxd = 0
	mxa = None
	
	n = 0

	for a1 in c1.cell:
		a2 = corresponding_atom(a1, c2)
		r = a1.distance(a2)
		d = r**2
		sd += d
		if mxd < d:
			mxd = d
			mxa = a1
		if a1.name() != 'X':
			n += 1
	print(n)
	sd = (sd / n) ** 0.5
	mxd = mxd ** 0.5
	return sd, mxd, mxa














