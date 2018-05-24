import re, os, sys, math
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys

lattice_str = 'DIRECT LATTICE VECTOR COMPONENTS (BOHR)'
basis_str = 'LOCAL ATOMIC FUNCTIONS BASIS SET'

atom_regex = '\s+([0-9]{1,3})\s+([A-Z][a-z]?)\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})'
orb_regex = '\s+(([0-9]+)\-\s+)?([0-9]+)\s+[SPDFGH]'

spin_regex = '\s+(ALPHA[\+\-]BETA)\s+ELECTRONS'
matrix_regex =  '\s+([A-Z]+\s+MATRIX)\s+\-\s+CELL\s+N\.\s+[0-9]+\(\s*(\-?[0-9]+)\s+(\-?[0-9]+)\s+(\-?[0-9]+)\s*\)'

NoSpin = 0
AplusB = 1
AminusB = 2

class CrystalMatrix:

	def __init__(self):
		self.cell = None

	@staticmethod
	def read_lattice(lines):
		for n in xrange(len(lines)):
			if lattice_str in lines[n]:
				break
		vectors = []
		for line in lines[n+1: n+4]:
			ls = line.split()
			assert len(ls) == 3
			v = Vector([(float(x) * Units.BOHR / Units.UNIT) for x in ls])
			vectors.append(v)
		return vectors

	@staticmethod
	def read_cell_and_basis(lines):
		atoms = []
		first_n = -1
		last_n = -1
		numorbmap = {}
		for n in xrange(len(lines)):
			if basis_str in lines[n]:
				break
		for line in lines[n+4:]:
			m = re.match(atom_regex, line)
			m1 = re.match(orb_regex, line)
			if not line:
				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
					al.data()[AtomKeys.ORBITAL_COUNT] = numorbmap[al.name()]
				break
			elif m:
				if atoms:
					al = atoms[-1]
					if last_n != -1 and first_n != -1:
						numorb = last_n - first_n + 1
						numorbmap[al.name()] = numorb
					al.data()[AtomKeys.ORBITAL_COUNT] = numorbmap[al.name()]
						
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
		return atoms


	@staticmethod
	def n_to_atomorb_tuple(k, atoms):
		for ca in atoms:
			no = ca.data()[AtomKeys.ORBITAL_COUNT]	
			if k <= no:
				return (ca.tuple_data(), k)
			else:
				k -= no
				continue
		assert False
			

	@staticmethod
	def read_matrix_small_part(lines, n, atoms, satoms, prec):
		data = {}
		ls = [int(i) for i in lines[n].split()]
		atomsv = atoms
		atomsh = satoms		#TODO: check if not and complex possible
		h_indices = [CrystalMatrix.n_to_atomorb_tuple(i, atomsh) for i in ls]
		ntot = reduce((lambda x,y: x + y), [ca.data()[AtomKeys.ORBITAL_COUNT] for ca in atoms])
		v_indices = [CrystalMatrix.n_to_atomorb_tuple(i+1, atomsv) for i in xrange(ntot)]
		n += 1
		if lines[n] == '':
			n += 1
		for k in xrange(n, n + ntot):
			ls = lines[k].split()
			assert len(ls) == len(h_indices) + 1
			vi = k - n + 1
			assert int(ls[0]) == vi, vi
			for vj in range(1, len(ls)):
				value = float(ls[vj])
				key = tuple(sorted([v_indices[vi - 1], v_indices[vj - 1]]))
				data[key] = value
		return data, k


	@staticmethod
	def read_matrix_cell_part(lines, n, cell, prec):
		print lines[n]
		m = re.match(matrix_regex, lines[n])
		assert m
		shifts = [m.group(i) for i in (2,3,4)]
		atoms = cell.cell
		satoms = cell.shift(shifts)
		n += 2
		datamap = {}
		while True:
			
			data, n = CrystalMatrix.read_matrix_small_part(lines, n, atoms, satoms, prec)
			datamap.update(data)
			n += 1
			if lines[n] != '':
				return datamap, n
			else:
				n += 1
				if re.match(matrix_regex, lines[n]) or re.match(spin_regex, lines[n]):
					return datamap, n
		
		

	@staticmethod
	def read_matrix(lines, n, cell, prec):
		datamap = {}
		while True:
			data, n = CrystalMatrix.read_matrix_cell_part(lines, n, cell, prec)
			if re.match(spin_regex, lines[n]):
				return datamap, n
			elif 'TTTTTTTTTTTTTT' in lines[n]:
				return datamap, n
		
			

	@staticmethod
	def read_matrices(lines, cell, prec):
		n = 0
		while True:
			mtype = NoSpin
			for k in xrange(n, len(lines)):
				m = re.match(spin_regex, lines[k])
				m1 = re.match(matrix_regex, lines[k])
				if m:
					mt = m.group(1)
					print lines[k]
					if mt == 'ALPHA+BETA':
						mtype = AplusB
					elif mt == 'ALPHA-BETA':
						mtype = AminusB
					else:
						assert False, mt
					k1 = k + 1
					for k in xrange(k1, len(lines)):
						m1 = re.match(matrix_regex, lines[k])
						if m1:
							break
					break
				elif m1:
					mtype = NoSpin
					break
			data, n = CrystalMatrix.read_matrix(lines, k, cell, prec)
			if 'TTTTTTTTTTTTTT' in lines[n]:
				break

			
	@staticmethod
	def from_string(datastring, prec):
		lines = datastring.splitlines()
		vectors = CrystalMatrix.read_lattice(lines)
		atoms = CrystalMatrix.read_cell_and_basis(lines)
		cm = CrystalMatrix()
		cm.cell = Cell(atoms, vectors)
		CrystalMatrix.read_matrices(lines, cm.cell, prec)
		return cm

	@staticmethod
	def from_file(name, prec=0):
		with open(name) as f:
			return CrystalMatrix.from_string(f.read(), prec)


