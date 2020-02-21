#!/usr/bin/python 
import re, os, sys, threading, time
from math3d import Vector
import fortranformat as ff
from ..util.units import Units
from ..structures.atom_vector import AtomVector
#TODO: merge with turbo_format if possible

GRAD_HEADER = '\s+cycle\s.*\sSCF energy\s.*\s\|dE/dxyz\|\s.*'

def read_atoms(filename):
	with open(filename) as inp:
		atoms = []
		lines = inp.read().splitlines()
		assert lines[0] == "$coord"
		assert lines[-1] == "$end"
		for line in lines[1:-1]:
			ls = line.strip().split()
			name = ls[3] if len(ls) == 4 else 'zz'
			vector = Vector([float(i) * Units.BOHR / Units.UNIT for i in ls[0:3]])
			atoms.append(AtomVector(name, vector))
		return atoms
		
def write_atoms(atoms, filename):
	with open(filename, 'w') as outp:
		outp.write("$coord\n")
		for atom in atoms:
			outp.write(" %16.10f %16.10f% 16.10f  %s\n" % (atom.position().x / Units.BOHR * Units.UNIT, atom.position().y / Units.BOHR * Units.UNIT, atom.position().z / Units.BOHR * Units.UNIT, atom.name()))
		outp.write("$end")

def get_last_grad(lines):
	res = 0
	for n in range(len(lines)):
		if re.match(GRAD_HEADER, lines[n]):
			res = n
	return res

def read_grad(atomnum):
	result = []
	with open("control") as c:
		lines = c.read().splitlines()
		n = get_last_grad(lines)
		for k in range(n + atomnum + 1, n + atomnum*2 + 1):
			result.extend([float(x.replace("D", "E")) for x in lines[k].split()])
		return result

def write_grad(x):
	newdata = ''
	atomnum = len(x) / 3
	with open("control") as c:
		lines = c.read().splitlines()
		n = get_last_grad(lines)
		for k in range(0, n + atomnum + 1):
			newdata += lines[k] + "\n"
		for k in range(n + atomnum + 1, n + atomnum*2 + 1):
			an = k - n - atomnum - 1
			newdata += ff.FortranRecordWriter("( 2X,D20.14,2X,D20.14,2X,D20.14 )").write(x[an*3:an*3+3]) + '\n'
			pass
		for k in range(n + atomnum*2 + 1, len(lines)):
			newdata += lines[k] + "\n"
	with open("control", 'w') as c:
		c.write(newdata)


def clear_relax():
	newdata = ''
	with open("control") as c:
		lines = c.read().splitlines()
		in_relax = False
		for line in lines:
			if '$forceapprox  ' in line:
				in_relax = True
			elif in_relax and '$' in line:
				in_relax = False	
			if not in_relax:
				newdata += line + "\n"
	with open("control", 'w') as c:
		c.write(newdata)
















