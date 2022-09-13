#!/usr/bin/python
import math, re
from math3d import Vector
from .units import Units
from .elements import ELEMENTS
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys

def atom_block(atoms, units=Units.ANGSTROM):
	res = ''
	k = Units.UNIT / units
	for a in atoms:
		res += ("%2s %17.10f %17.10f %17.10f\n" % (a.name(), a.position().x * k, a.position().y * k, a.position().z * k))
	return res

def write_xyz(atoms, name):
	res = atom_block(atoms)
	with open(name, 'w') as f:
		f.write("%d\n\n" % len(atoms))
		f.write(res)

def read_xyz(name):
	nat = -1
	with open(name, 'r') as f:
		lines = f.read().splitlines()
		n = 0
		atoms = []
		if re.match('^[0-9]+$', lines[0].strip()):
			nat = int(lines[0].strip())
			n = 2
		while n < len(lines):
			ls = lines[n].split()
			if len(ls) != 4 or nat == 0:
				break
			name = ls[0]
			if re.match('^[0-9]+$', name):
				num = int(name)
				name = ELEMENTS[num].symbol
			v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in ls[1:]])
			a = AtomVector(name, v)
			atoms.append(a)
			n += 1
			nat -= 1
		cell = Cell(atoms, [], None, None, None)
		return cell
			
		












