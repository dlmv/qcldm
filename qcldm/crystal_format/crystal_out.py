import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
from ..functions.gauss_function import GaussFunctionContracted, GaussFunctionNormed

DFT_PARAMS = 'DFT PARAMETERS'

BASIS = 'LOCAL ATOMIC FUNCTIONS BASIS SET'

dft_regex = '\s+[0-9]+\s+([0-9]+)\s+([A-Z][a-z]?)\s+([0-9]{1,3}\.[0-9]{4})\s+[0-9\.\-]+'

atom_regex = '\s+[0-9]+\s+[0-9]+\s+([A-Z][a-z]?)' + '\s+(\-?[0-9]+\.[0-9]+E[\+\-][0-9]+)' * 3
#                 #            n         type          xyz      


basis_atom_regex = '\s+([0-9]{1,3})\s+([A-Z][a-z]?)\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})'
basis_orb_regex = '\s+(([0-9]+)\-\s+)?([0-9]+)\s+([SPDFGH])'
basis_gauss_regex = '\s{20}' + '\s+(\-?[0-9]+\.[0-9]+E[\+\-][0-9]+)' * 4

ATOMS_POP = 'ATOM    Z CHARGE  SHELL POPULATION'

pop_regex = '\s+([0-9]+)\s+([A-Z][a-z]?)\s+[0-9]+\s+([0-9]{1,3}\.[0-9]{3})\s+'
#                 n         type           z            charge

GEOM_OUT = 'GEOMETRY OUTPUT FILE'
LATTICE = 'DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)'
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
			for cg in self.basis[a.name()]:
				cut = max(cut, cg.get_cutoff(prec))
			a.data()[AtomKeys.CUTOFF] = cut

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
	def get_cell(lines, vm, basis):
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
				numorb = 0
				for cg in basis[a.name()]:
					numorb += cg.fs[0][1].l * 2 + 1
				a.data()[AtomKeys.ORBITAL_COUNT] = numorb
				a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
				a.data()[AtomKeys.ESTIMATED_VALENCE] = Shells.estimate_valence_byname(a.name())
				
				atoms.append(a)
		cell = Cell(atoms, vectors)
		return cell
		
	@staticmethod
	def read_basis(lines):
		numorbmap = {}
		basismap = {}
		basis = []
		atoms = []
		first_n = -1
		last_n = -1
		l = 0
		for n in xrange(len(lines)):
			if BASIS in lines[n]:
				break
				
		gs = []
		for line in lines[n+4:]:
			m = re.match(basis_atom_regex, line)
			m1 = re.match(basis_orb_regex, line)
			m2 = re.match(basis_gauss_regex, line)
			if not line:
				if gs:
					cg = GaussFunctionContracted()
					cg.fs = gs
					basis.append(cg)
				gs = []

				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
					if basis:
						print basis
						basismap[al.name()] = basis
						basis = []
				break
			elif m:
				if gs:
					cg = GaussFunctionContracted()
					cg.fs = gs
					basis.append(cg)
				gs = []
				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
					if basis:
						basismap[al.name()] = basis
						basis = []
						
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
				lname = m1.group(4)
				l = Shells.SHELLS.index(lname)
				if gs:
					cg = GaussFunctionContracted()
					cg.fs = gs
					basis.append(cg)
				gs = []
			elif m2:
				ls = line.split()
				a = float(ls[0])
				ln = min(l + 1, 3)
				k = float(ls[ln])
				g = GaussFunctionNormed(a, l)
				gs.append([k, g])
				assert k
		
#		for k in basismap.keys():
#			s = 0
#			for cg in basismap[k]:
#				s += cg.fs[0][1].l * 2 + 1
#			print k, s, numorbmap[k]
			
				
		return basismap
			

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


