#!/usr/bin/python

import sys, os, logging

from ..openmx_format.dat_format import DAT_INPUT
from cell_neighbours import distance, angle

def corresponding_atom(a, cell2):
	return cell2.cell[a.num - 1].shifted(a.shifts)

def compare_cells(c1, c2):
	if len(c1.cell) != len(c2.cell):
		print 'different atom number!!'
		return
	for a1, a2 in zip(c1.cell, c2.cell):
		if a1.name() != a2.name():
			print 'different atom types!'
			return
	c1.neighbours.load_all()
	c2.neighbours.load_all()
	c1.neighbours.merge_from(c2.neighbours)

	sd = 0
	mxd = 0

	n = 0

	a, nb = 0,0

	for a1 in c1.cell:
		a2 = corresponding_atom(a1, c2)
		for n1 in c1.neighbours.first_neighbours(a1):
			n2 = corresponding_atom(n1, c2)
			r1 = distance(a1, n1)
			r2 = distance(a2, n2)
			d = ((r1 - r2) / r1)**2
			sd += d
			if mxd < d and n1 in c1.cell:
				mxd = d
				a, nb = a1, n1
			n += 1
	sd = (sd / n) ** 0.5
	mxd = mxd ** 0.5
	return sd, mxd, a, nb
















