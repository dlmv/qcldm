import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
from ..functions.gauss_function import GaussFunctionContracted, GaussFunctionNormed

DFT_PARAMS = 'DFT PARAMETERS'

BASIS = 'LOCAL ATOMIC FUNCTIONS BASIS SET'

dft_regex = '\s+[0-9]+\s+([0-9]+)\s+([A-Z][a-zA-Z]?)\s+([0-9]{1,3}\.[0-9]{4})\s+(\-?[0-9]{1,3}\.[0-9]{4})\s+[0-9\.\-]+'

atom_regex = '\s+[0-9]+\s+[0-9]+\s+([A-Z][a-zA-Z]?)' + '\s+(\-?[0-9]+\.[0-9]+E[\+\-][0-9]+)' * 3
#                 #            n         type          xyz      
atom1_regex = '\s+[0-9]+\s+T\s+[0-9]+\s+([A-Z][a-zA-Z]?)' + '\s+(\-?[0-9]+\.[0-9]+E[\+\-][0-9]+)' * 3

basis_atom_regex = '\s+([0-9]{1,3})\s+([A-Z][a-zA-Z]?)\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})'
basis_orb_regex = '\s+(([0-9]+)\-\s+)?([0-9]+)\s+([SPDFGH]{1,2})'
basis_gauss_regex = '\s{20}' + '\s*([\- ][0-9]+\.[0-9]+E[\+\-][0-9]+)' * 4

ATOMS_POP = 'ATOM    Z CHARGE  SHELL POPULATION'

pop_regex = '\s+([0-9]+)\s+([A-Z][a-zA-Z]?)\s+[0-9]+\s+([0-9]{1,3}\.[0-9]{3})\s+'
#                 n         type           z            charge


#float_regex = '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'


GEOM_OUT = 'GEOMETRY OUTPUT FILE'
LATTICE = 'DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)'
NOLATTICE = '(NON PERIODIC DIRECTION: LATTICE PARAMETER FORMALLY SET TO'
COORDS =  'CARTESIAN COORDINATES - PRIMITIVE CELL'


class CrystalOut:

	def __init__(self):
		self.name = None
		self.cell = None
		self.valence_map = {}
		self.basis = {}
	
	def get_cutoffs(self, prec=0.0001):
		cut = 0
		for a in self.cell.atoms:
			for l in self.basis[a.name()].keys():
				for cg in self.basis[a.name()][l]:
					cut = max(cut, cg.get_cutoff(prec))
			a.data()[AtomKeys.CUTOFF] = cut * Units.BOHR / Units.UNIT

	@staticmethod
	def get_charge_map(lines):
		vm = {}
		for n in xrange(len(lines)):
			if DFT_PARAMS in lines[n]:
				break
		n += 3
		for k in xrange(n, len(lines)):
			if not lines[k]:
				break
			m = re.match(dft_regex, lines[k])
			if m:
				vm[CrystalOut.get_normal_name(m.group(2))] = int(float(m.group(3)) + float(m.group(4)))
		return vm

	@staticmethod
	def get_cell(lines, vm, basis):
		k = 0
		for n in xrange(len(lines)):
			if GEOM_OUT in lines[n]:
				k = n
				break
		n = k
		for k in xrange(n, len(lines)):
			if NOLATTICE in lines[k]:
				break
		if k == len(lines) - 1:
			k = 0
		for n in xrange(k, len(lines)):
			if LATTICE in lines[n]:
				k = n
				break
		vectors = []
		atoms = []
		if LATTICE in lines[k]:
			for n in xrange(k+2, k+5):
				ls = lines[n].split()
				assert len(ls) == 3
				v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in ls])
				vectors.append(v)

			for k in xrange(n, len(lines)):
				if COORDS in lines[k]:
					break
			for n in xrange(k+4, len(lines)):
				m = re.match(atom_regex, lines[n])
			
				if m:
					v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in [m.group(2), m.group(3), m.group(4)]])
					name = CrystalOut.get_normal_name(m.group(1))
					a = AtomVector(name, v)
					numorb = 0
					for l in basis[a.name()].keys():
						for cg in basis[a.name()][l]:
							numorb += cg.fs[0][1].l * 2 + 1
					a.data()[AtomKeys.ORBITAL_COUNT] = numorb
					a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
					atoms.append(a)
		elif NOLATTICE in lines[k]:
#			vectors = [Vector([500 if i == j else 0 for j in range(3)]) for i in range(3)]
			for n in xrange(k+5, len(lines)):
				m = re.match(atom1_regex, lines[n])
			
				if m:
					v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in [m.group(2), m.group(3), m.group(4)]])
					name = CrystalOut.get_normal_name(m.group(1))
					a = AtomVector(name, v)
					numorb = 0
					for l in basis[a.name()].keys():
						for cg in basis[a.name()][l]:
							numorb += cg.fs[0][1].l * 2 + 1
					a.data()[AtomKeys.ORBITAL_COUNT] = numorb
					a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
					atoms.append(a)


		cell = Cell(atoms, vectors)
		return cell

#	@staticmethod
#	def read_input_basis(lines):
#		basismap = {}
#		basis = {}
#		for n in xrange(len(lines)):
#			if 'CRYSTAL' in lines[n]:
#				break
		
	@staticmethod
	def read_basis(lines):
		numorbmap = {}
		basismap = {}
		basis = {}
		atoms = []
		first_n = -1
		last_n = -1
		l = 0
		l1 = -1
		for n in xrange(len(lines)):
			if BASIS in lines[n]:
				break
				
		gs = []
		gs1 = []
		for line in lines[n+4:]:
			m = re.match(basis_atom_regex, line)
			m1 = re.match(basis_orb_regex, line)
			m2 = re.match(basis_gauss_regex, line)
			if not line:
				if gs:
					cg = GaussFunctionContracted()
					cg.fs = gs
					if gs[0][1].l not in basis.keys():
						basis[gs[0][1].l] = []
					basis[gs[0][1].l].append(cg)
				if gs1:
					cg = GaussFunctionContracted()
					cg.fs = gs1
					if gs1[0][1].l not in basis.keys():
						basis[gs1[0][1].l] = []
					basis[gs1[0][1].l].append(cg)
				gs = []
				gs1 = []

				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
					if basis:
						basismap[al.name()] = basis
						basis = {}
				break
			elif m:
				if gs:
					cg = GaussFunctionContracted()
					cg.fs = gs
					if gs[0][1].l not in basis.keys():
						basis[gs[0][1].l] = []
					basis[gs[0][1].l].append(cg)
				if gs1:
					cg = GaussFunctionContracted()
					cg.fs = gs1
					if gs1[0][1].l not in basis.keys():
						basis[gs1[0][1].l] = []
					basis[gs1[0][1].l].append(cg)
				gs = []
				gs1 = []
				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
					if basis:
						basismap[al.name()] = basis
						basis = {}
						
				v = Vector([(float(x) * Units.BOHR / Units.UNIT) for x in [m.group(3), m.group(4), m.group(5)]])
				name = CrystalOut.get_normal_name(m.group(2))
				a = AtomVector(name, v)
				atoms.append(a)
				first_n = -1
				last_n = -1
			elif m1:
				if first_n == -1:
					if m1.group(1):
						first_n = int(m1.group(2))
					else:
						first_n = int(m1.group(3))
				last_n = int(m1.group(3))
				lname = m1.group(4)
				if lname == 'SP':
					l = 0
					l1 = 1
				elif len(lname) == 1:
					l = Shells.SHELLS.index(lname)
					l1 = -1
				else:
					assert False
				if gs:
					cg = GaussFunctionContracted()
					cg.fs = gs
					if gs[0][1].l not in basis.keys():
						basis[gs[0][1].l] = []
					basis[gs[0][1].l].append(cg)
				if gs1:
					cg = GaussFunctionContracted()
					cg.fs = gs1
					if gs1[0][1].l not in basis.keys():
						basis[gs1[0][1].l] = []
					basis[gs1[0][1].l].append(cg)
				gs = []
				gs1 = []
			elif m2:
				ls = line.split()
				a = float(m2.group(1))
				ln = min(l + 1, 3)
				k = float(m2.group(ln + 1))
				g = GaussFunctionNormed(a, l)
				gs.append([k, g])
				if l1 != -1:
					ln = min(l1 + 1, 3)
					k = float(m2.group(ln + 1))
					g = GaussFunctionNormed(a, l1)
					gs1.append([k, g])	
				assert k
		
#		for k in basismap.keys():
#			s = 0
#			for cg in basismap[k]:
#				s += cg.fs[0][1].l * 2 + 1
#			print k, s, numorbmap[k]


		return basismap
		
	@staticmethod
	def get_normal_name(a):
		return a[0].upper() + a[1:].lower()
			

	@staticmethod
	def get_mulliken_pop(lines, cell):
		pop = {}
		for n in xrange(len(lines)):
			if ATOMS_POP in lines[n]:
				break
		n += 1
		while not lines[n]:
			n += 1
		for k in xrange(n, len(lines)):
			if not lines[k]:
				break
			m = re.match(pop_regex, lines[k])
			
			if m:
				pop[CrystalOut.get_normal_name(m.group(2))] = float(m.group(3))
		for a in cell.atoms:
			a.data()[AtomKeys.MULLIKEN_CHARGE] = a.data()[AtomKeys.FULL_VALENCE] - pop[a.name()]
		
		
	@staticmethod
	def from_string(datastring, name):
		lines = datastring.splitlines()
		co = CrystalOut()
		co.name = name
		co.valence_map = CrystalOut.get_charge_map(lines)
		co.basis = CrystalOut.read_basis(lines)
		co.cell = CrystalOut.get_cell(lines, co.valence_map, co.basis)
		CrystalOut.get_mulliken_pop(lines, co.cell)
		logging.debug(u'ATOM VALENCE:')
		for k in co.valence_map.keys():
			logging.debug(u'%s: %d' % (k, co.valence_map[k]))
		return co
			

	@staticmethod
	def from_file(name):
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Reading output from %s' % name)
		logging.info(u'*********************************************')
		logging.info(u'')
		with open(name) as f:
			return CrystalOut.from_string(f.read(), name)


