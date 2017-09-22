#!/usr/bin/python
import math
from units import Units

def atom_block(atoms, units=Units.ANGSTROM):
	res = ''
	k = Units.UNIT / units
	for a in atoms:
		res += ("%2s %15.8f %15.8f %15.8f\n" % (a.name(), a.position().x * k, a.position().y * k, a.position().z * k))
	return res

def write_xyz(atoms, name):
	res = atom_block(atoms)
	with open(name, 'w') as f:
		f.write("%d\n\n" % len(atoms))
		f.write(res)














