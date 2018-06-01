import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells

DFT_PARAMS = 'DFT PARAMETERS'

BASIS = 'LOCAL ATOMIC FUNCTIONS BASIS SET'

dft_regex = '\s+[0-9]+\s+([0-9]+)\s+([A-Z][a-z]?)\s+([0-9]{1,3}\.[0-9]{4})\s+[0-9\.\-]+'

atom_regex = '\s+[0-9]+\s+[0-9]+\s+([A-Z][a-z]?)' + '\s+(\-?[0-9]+\.[0-9]+E[\+\-][0-9]+)' * 3
#                 #            n         type          xyz      


basis_atom_regex = '\s+([0-9]{1,3})\s+([A-Z][a-z]?)\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})'
basis_orb_regex = '\s+(([0-9]+)\-\s+)?([0-9]+)\s+[SPDFGH]'

ATOMS_POP = 'ATOM    Z CHARGE  SHELL POPULATION'

pop_regex = '\s+([0-9]+)\s+([A-Z][a-z]?)\s+[0-9]+\s+([0-9]{1,3}\.[0-9]{3})\s+'
#                 n         type           z            charge

GEOM_OUT = 'GEOMETRY OUTPUT FILE'
LATTICE = 'DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)'
COORDS =  'CARTESIAN COORDINATES - PRIMITIVE CELL'


class CrystalOut:

	def __init__(self):
		self.cell = None
		self.valence_map = {}

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
				vm[m.group(2)] = int(float(m.group(3)))
		return vm

	@staticmethod
	def get_cell(lines, vm, nm):
		for n in xrange(len(lines)):
			if GEOM_OUT in lines[n]:
				break
		for k in xrange(n, len(lines)):
			if LATTICE in lines[k]:
				break
		vectors = []
		for n in xrange(k+2, k+5):
			ls = lines[n].split()
			assert len(ls) == 3
			v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in ls])
			vectors.append(v)

		for k in xrange(n, len(lines)):
			if COORDS in lines[k]:
				break
		atoms = []
		for n in xrange(k+4, len(lines)):
			m = re.match(atom_regex, lines[n])
			
			if m:
				v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in [m.group(2), m.group(3), m.group(4)]])
				name = m.group(1)
				a = AtomVector(name, v)
				a.data()[AtomKeys.ORBITAL_COUNT] = nm[a.name()]
				a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
				a.data()[AtomKeys.ESTIMATED_VALENCE] = Shells.estimate_valence_byname(a.name())
				
				atoms.append(a)
		cell = Cell(atoms, vectors)
		return cell
		
	@staticmethod
	def read_orbnums(lines):
		numorbmap = {}
		atoms = []
		first_n = -1
		last_n = -1
		for n in xrange(len(lines)):
			if BASIS in lines[n]:
				break
		for line in lines[n+4:]:
			m = re.match(basis_atom_regex, line)
			m1 = re.match(basis_orb_regex, line)
			if not line:
				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
				break
			elif m:
				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
						
				v = Vector([(float(x) * Units.BOHR / Units.UNIT) for x in [m.group(3), m.group(4), m.group(5)]])
				name = m.group(2)
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
		return numorbmap
			

	@staticmethod
	def get_mulliken_pop(lines, cell):
		pop = {}
		for n in xrange(len(lines)):
			if ATOMS_POP in lines[n]:
				break
		n += 3
		for k in xrange(n, len(lines)):
			if not lines[k]:
				break
			m = re.match(pop_regex, lines[k])
			if m:
				pop[m.group(2)] = float(m.group(3))
		for a in cell.atoms:
			a.data()[AtomKeys.MULLIKEN_CHARGE] = a.data()[AtomKeys.FULL_VALENCE] - pop[a.name()]
		
		
	@staticmethod
	def from_string(datastring):
		lines = datastring.splitlines()
		co = CrystalOut()

		co.valence_map = CrystalOut.get_charge_map(lines)
		nm = CrystalOut.read_orbnums(lines)
		co.cell = CrystalOut.get_cell(lines, co.valence_map, nm)
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
			return CrystalOut.from_string(f.read())


