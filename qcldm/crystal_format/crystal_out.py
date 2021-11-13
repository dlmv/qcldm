import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
from ..gauss_functions.gauss_function import GaussFunctionContracted, GaussFunctionNormed

DFT_PARAMS = 'DFT PARAMETERS'

BASIS = 'LOCAL ATOMIC FUNCTIONS BASIS SET'

dft_regex = '\s+[0-9]+\s+([0-9]+)\s+([A-Z][a-zA-Z]?)\s+([0-9]{1,3}\.[0-9]{4})\s+(\-?[0-9]{1,3}\.[0-9]{4})\s+[0-9\.\-]+'

atom_regex = '\s+[0-9]+\s+[0-9]+\s+([A-Z][a-zA-Z]?)' + '\s+(\-?[0-9]+\.[0-9]+[Ee][\+\-][0-9]+)' * 3
#                 #            n         type          xyz      
atom1_regex = '\s+[0-9]+\s+T\s+[0-9]+\s+([A-Z][a-zA-Z]?)' + '\s+(\-?[0-9]+\.[0-9]+[Ee][\+\-][0-9]+)' * 3

basis_atom_regex = '\s+([0-9]{1,3})\s+([A-Z][a-zA-Z]?)\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})'
basis_orb_regex = '\s+(([0-9]+)\-\s+)?([0-9]+)\s+([SPDFGH]{1,2})'
basis_gauss_regex = '\s{20}' + '\s*([\- ][0-9]+\.[0-9]+[Ee][\+\-][0-9]+)' * 4

ATOMS_POP = 'ATOM    Z CHARGE  SHELL POPULATION'

pop_regex = '\s+([0-9]+)\s+([A-Z][a-zA-Z]?)\s+[0-9]+\s+([0-9]{1,3}\.[0-9]{3})\s+'
#                 n         type           z            charge


#float_regex = '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'


GEOM_OUT = 'GEOMETRY OUTPUT FILE'
LATTICE_MAT = 'DIRECT LATTICE VECTORS CARTESIAN COMPONENTS (ANGSTROM)'
LATTICE_DEG = 'LATTICE PARAMETERS (ANGSTROMS AND DEGREES)'
NOLATTICE = '(NON PERIODIC DIRECTION: LATTICE PARAMETER FORMALLY SET TO'
ASSYM = ' ATOMS IN THE ASYMMETRIC UNIT '
COORDS =  'CARTESIAN COORDINATES - PRIMITIVE CELL'
TRANS_MATRIX = 'TRANSFORMATION MATRIX PRIMITIVE-CRYSTALLOGRAPHIC CELL'


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
	def get_geom_start(lines):
		for n in xrange(len(lines)):
			if GEOM_OUT in lines[n]:
				return n
		return 0

	@staticmethod
	def get_last_geom(lines):
		n = 0
		for k in xrange(len(lines)):
			if NOLATTICE in lines[k]:
				n = k
		return n

	@staticmethod
	def get_first_geom(lines):
		for k in xrange(len(lines)):
			if NOLATTICE in lines[k]:
				return k
		return 0

	
		

	@staticmethod
	def get_cell_new(lines, vm, basis):
		start = CrystalOut.get_first_geom(lines)
		if ASSYM in lines[start + 2]:
			atoms = []
			logging.debug(u'NON-PERIODIC SYSTEM')
			for n in xrange(start+5, len(lines)):
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
			
			return Cell(atoms, [], None, [], range(len(atoms)))
		elif LATTICE_DEG in lines[start + 2]:
			atoms = []
			vectors = []
			assym_n = []
			trans_mat = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
			logging.debug(u'PERIODIC SYSTEM')
			m = re.match('\s*ATOMS IN THE ASYMMETRIC UNIT\s+([0-9]+)\s+- ATOMS IN THE UNIT CELL:\s+([0-9]+)', lines[start+7])
			ass, tot = [int(x) for x in m.groups()]
			for n in range(tot):
				if ' T ' in lines[start + 10 + n]:
					assym_n.append(n)
			symstart = 0
			if TRANS_MATRIX in lines[start + tot + 11]:
				logging.debug(u'PRIMITIVE CELL IS DIFFERENT FROM CRYSTALLOGRAPHIC')
				mat_raw = [float(i) for i in lines[start + tot + 12].split()]
				trans_mat = map(lambda x: mat_raw[3*x:(x+1)*3], range(3))
				symstart = start + 2*tot + 25
			else:
				logging.debug(u'PRIMITIVE CELL IS EQUAL TO CRYSTALLOGRAPHIC')
				symstart = start + tot + 13
			m = re.match('\s+\*+\s+([0-9]+)\s+SYMMOPS - TRANSLATORS IN FRACTIONAL UNITS', lines[symstart])
			nops = int(m.group(1))
			symops = []
			for i in range(nops):
				ls = lines[symstart + 3 + i].split()
				n, inv = int(ls[0]), int(ls[1])
				rot_raw = [float(i) for i in ls[2:11]]
				rot_mat = map(lambda x: rot_raw[3*x:(x+1)*3], range(3))
				translator = [float(i) for i in ls[11:14]]
				symop = [n,inv,rot_mat,translator]
				symops.append(symop)
			for n in range(symstart + nops + 6, symstart + nops + 9):
				ls = lines[n].split()
				assert len(ls) == 3
				v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in ls])
				vectors.append(v)
			coordstart = symstart + nops + 15
			for n in range(coordstart, coordstart + tot):
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
					if vm:
						a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
					atoms.append(a)
			return Cell(atoms, vectors, trans_mat, symops, assym_n)
		
	
	
		

#	@staticmethod
#	def get_cell(lines, vm, basis):
#		k = 0
#		for n in xrange(len(lines)):
#			if GEOM_OUT in lines[n]:
#				k = n
#				break
#		n = k
#		for k in xrange(n, len(lines)):
#			if NOLATTICE in lines[k]:
#				break
#		if k == len(lines) - 1:
#			k = 0
#		for n in xrange(k, len(lines)):
#			if LATTICE_MAT in lines[n]:
#				k = n
#				break
#		vectors = []
#		atoms = []
#		if LATTICE_MAT in lines[k]:
#			for n in xrange(k+2, k+5):
#				ls = lines[n].split()
#				assert len(ls) == 3
#				v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in ls])
#				vectors.append(v)

#			for k in xrange(n, len(lines)):
#				if COORDS in lines[k]:
#					break
#			for n in xrange(k+4, len(lines)):
#				m = re.match(atom_regex, lines[n])
#			
#				if m:
#					v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in [m.group(2), m.group(3), m.group(4)]])
#					name = CrystalOut.get_normal_name(m.group(1))
#					a = AtomVector(name, v)
#					numorb = 0
#					for l in basis[a.name()].keys():
#						for cg in basis[a.name()][l]:
#							numorb += cg.fs[0][1].l * 2 + 1
#					a.data()[AtomKeys.ORBITAL_COUNT] = numorb
#					a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
#					atoms.append(a)
#		elif NOLATTICE in lines[k]:
##			vectors = [Vector([500 if i == j else 0 for j in range(3)]) for i in range(3)]
#			for n in xrange(k+5, len(lines)):
#				m = re.match(atom1_regex, lines[n])
#			
#				if m:
#					v = Vector([(float(x) * Units.ANGSTROM / Units.UNIT) for x in [m.group(2), m.group(3), m.group(4)]])
#					name = CrystalOut.get_normal_name(m.group(1))
#					a = AtomVector(name, v)
#					numorb = 0
#					for l in basis[a.name()].keys():
#						for cg in basis[a.name()][l]:
#							numorb += cg.fs[0][1].l * 2 + 1
#					a.data()[AtomKeys.ORBITAL_COUNT] = numorb
#					a.data()[AtomKeys.FULL_VALENCE] = vm[a.name()]
#					atoms.append(a)


#		cell = Cell(atoms, vectors)
#		return cell

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
		for n in xrange(len(lines)):
			if ATOMS_POP in lines[n]:
				break
		n += 1
		if n >= len(lines):
			return
		while not lines[n]:
			n += 1
		for k in xrange(n, len(lines)):
			if not lines[k]:
				break
			m = re.match(pop_regex, lines[k])
			
			if m:
				atomnum = int(m.group(1)) - 1
				a = cell.atoms[atomnum]
				a.data()[AtomKeys.MULLIKEN_CHARGE] = a.data()[AtomKeys.FULL_VALENCE] - float(m.group(3))
		
		
	@staticmethod
	def from_string(datastring, name):
		lines = datastring.splitlines()
		co = CrystalOut()
		co.name = name
		co.valence_map = CrystalOut.get_charge_map(lines)
		co.basis = CrystalOut.read_basis(lines)
		co.cell = CrystalOut.get_cell_new(lines, co.valence_map, co.basis)
#		co.cell = CrystalOut.get_cell(lines, co.valence_map, co.basis)
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


