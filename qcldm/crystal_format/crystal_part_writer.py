import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..util.elements import ELEMENTS
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
from ..gauss_functions.gauss_function import GaussFunctionContracted, GaussFunctionNormed
import numpy as np
from ..atom.shells import Shells


DIVIDER = ' *******************************************************************************\n'
PARAMS_HEADER = '         A              B              C           ALPHA      BETA       GAMMA \n'

def unit_vector(vector):
	return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
	v1_u = unit_vector(v1)
	v2_u = unit_vector(v2)
	return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def angle_list(vectors):
	return [180. / math.pi * angle_between(v1, v2) for v1, v2 in zip([vectors[-1]] + list(vectors)[:-1], list(vectors)[1:] + [vectors[0]])]

def write_structure(cell, f):
	crysmat = np.array(cell.cryst_mat)
	primbas = np.array([v._data for v in cell.vectors])
	crysbas = np.transpose(crysmat).dot(primbas)
	f.write(' (NON PERIODIC DIRECTION: LATTICE PARAMETER FORMALLY SET TO 500)\n')
	f.write(DIVIDER)
	f.write(' LATTICE PARAMETERS (ANGSTROMS AND DEGREES) - BOHR = 0.5291772083 ANGSTROM\n')
	f.write(' PRIMITIVE CELL - CENTRING CODE 0/0 VOLUME= %12.6f - DENSITY %6.3f g/cm^3\n' % ((np.cross(primbas[0], primbas[1])).dot(primbas[2]), 0))
	f.write(PARAMS_HEADER)
	f.write(' %14.8f %14.8f %14.8f   %9.6f %9.6f %9.6f\n' % tuple([np.linalg.norm(v) for v in primbas] + angle_list(list(primbas))))
	f.write(DIVIDER)
	f.write(' ATOMS IN THE ASYMMETRIC UNIT    %d - ATOMS IN THE UNIT CELL:   %d\n' % (len(cell.assym_n), len(cell.atoms)))
	f.write('     ATOM                 X/A                 Y/B                 Z/C    \n')
	f.write(DIVIDER)	
	for i, atom in enumerate(cell.atoms):
		relatom = cell.relative_atom(atom, False)
		f.write('    %3d %s %3d %2s   %19.12e %19.12e %19.12e\n' % (i+1,
			'T' if i in cell.assym_n else 'F',
			ELEMENTS[relatom.name()].number,
			atom.name() + " " if len(atom.name()) == 1 else atom.name(),
			relatom.position()[0], relatom.position()[1], relatom.position()[2]
			))
	f.write("\n")
	if not (crysmat == np.identity(3)).all():
		f.write(" TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL\n")
		fs = ' %8.4f'* 9 + "\n"
		f.write(fs % tuple([v for row in cell.cryst_mat for v in row]))
		f.write("\n")
		f.write(DIVIDER)
		f.write(' CRYSTALLOGRAPHIC CELL (VOLUME=     %15.8f)\n' %  (np.cross(crysbas[0], crysbas[1])).dot(crysbas[2]))
		f.write(PARAMS_HEADER)
		f.write(' %14.8f %14.8f %14.8f   %9.6f %9.6f %9.6f\n' % tuple([np.linalg.norm(v) for v in crysbas] + angle_list(list(crysbas))))
		f.write("\n")
		f.write(' COORDINATES IN THE CRYSTALLOGRAPHIC CELL\n')
		f.write('     ATOM                 X/A                 Y/B                 Z/C    \n')
		f.write(DIVIDER)
		for i, atom in enumerate(cell.atoms):
			relatom = cell.relative_atom(atom, True)
			f.write('    %3d %s %3d %2s   %19.12e %19.12e %19.12e\n' % (i+1,
				'T' if i in cell.assym_n else 'F',
				ELEMENTS[relatom.name()].number,
				relatom.name() + " " if len(relatom.name()) == 1 else '',
				relatom.position()[0], relatom.position()[1], relatom.position()[2]
				))
		f.write("\n")
	f.write(" T = ATOM BELONGING TO THE ASYMMETRIC UNIT\n")
	f.write("\n")
	f.write(' ****  %3d SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS\n' % len(cell.symops))
	f.write(' **** MATRICES AND TRANSLATORS IN THE CRYSTALLOGRAPHIC REFERENCE FRAME\n')
	f.write('   V INV                    ROTATION MATRICES                   TRANSLATORS\n')
	for n,inv,rot_mat,translator in cell.symops:
		fs = ' %3d' * 2 + ' %5.2f' * 12 + "\n"
		f.write(fs % tuple([n, inv] + [v for row in rot_mat for v in row] + translator))
	f.write("\n")
	f.write(" DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)\n")
	f.write("          X                    Y                    Z\n")
	for v in cell.vectors:
		f.write("   %19.12e  %19.12e  %19.12e\n" % (v.x, v.y, v.z))
	f.write("\n")
	f.write("\n")
	f.write(' CARTESIAN COORDINATES - PRIMITIVE CELL\n')
	f.write(DIVIDER)
	f.write(" *      ATOM          X(ANGSTROM)         Y(ANGSTROM)         Z(ANGSTROM)\n")
	f.write(DIVIDER)

	for i, atom in enumerate(cell.atoms):
		f.write('  %3d %s %3d %2s   %19.12e %19.12e %19.12e\n' % (i+1,
			' ',
			ELEMENTS[atom.name()].number,
			atom.name() + " " if len(atom.name()) == 1 else atom.name(),
			atom.position().x, atom.position().y, atom.position().z
			))
	f.write("\n")
	f.write(DIVIDER)

def write_valence(co, f):
	f.write("\n")
	f.write(' DFT PARAMETERS\n')
	f.write("\n")
	f.write('     ATOM       ELECTRONS   NET CHARGE   R(ANGSTROM)\n')
	for name in list(co.valence_map.keys()):
		name2 = name + ' ' if len(name) == 1 else name
		f.write('   0 %3d  %s   %9.4f      0.0000     0.00000000\n' % (ELEMENTS[name].number, name2, co.valence_map[name]))
	f.write("\n")
	
def write_pop(cell, f):
	f.write("\n")
	f.write('  ATOM    Z CHARGE  SHELL POPULATION\n')
	f.write("\n")
	for i, atom in enumerate(cell.atoms):
		f.write(' %3d %s %3d %6.3f   \n' % (i+1,
			atom.name() + " " if len(atom.name()) == 1 else atom.name(),
			ELEMENTS[atom.name()].number,
			atom.data()[AtomKeys.FULL_VALENCE] - atom.data()[AtomKeys.MULLIKEN_CHARGE]
			))
	f.write("\n")

def write_basis(co, f):
	f.write("\n")
	f.write(DIVIDER)
	f.write(' LOCAL ATOMIC FUNCTIONS BASIS SET\n')
	f.write(DIVIDER)
	f.write('   ATOM   X(AU)   Y(AU)   Z(AU)  N. TYPE  EXPONENT  S COEF   P COEF   D/F/G COEF\n')
	f.write(DIVIDER)
	n = 1
	norb = {}
	for i, atom in enumerate(co.cell.atoms):

		name2 = atom.name() + " " if len(atom.name()) == 1 else atom.name()
		f.write(" %3d %2s  %6.3f  %6.3f  %6.3f\n" % (i+1, name2, atom.position().x / Units.BOHR, atom.position().y / Units.BOHR, atom.position().z / Units.BOHR))
		basis = co.basis[atom.name()]
		if atom.name() in list(norb.keys()):
			n += norb[atom.name()]
			continue
		k = n
		for l in sorted(basis.keys()):
			for cg in basis[l]:
				if l == 0:
					f.write('                               %4d S  \n' % (n))
				else:
					f.write('                        %4d-  %4d %s  \n' % (n,n + (2 * l),Shells.SHELLS[l]))
				n += (2 * l + 1)
				for c,g in cg.fs:
					ls = [g.a, 0., 0., 0.]
					ls[min(3, l + 1)] = c
					f.write('                                         %9.3e %9.3e %9.3e %9.3e\n' % tuple(ls))
		norb[atom.name()] = n - k
					

def write_crystal_part(co, fname):
	with open(fname, 'w') as f:
		write_structure(co.cell, f)
		write_basis(co, f)
		write_valence(co, f)
		write_pop(co.cell, f)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		

