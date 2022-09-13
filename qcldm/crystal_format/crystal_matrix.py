import re, os, sys, math, logging
from collections import Counter
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys
from ..atom.shells import Shells
from functools import reduce

#lattice_str = 'DIRECT LATTICE VECTOR COMPONENTS (BOHR)'
#basis_str = 'LOCAL ATOMIC FUNCTIONS BASIS SET'

#atom_regex = '\s+([0-9]{1,3})\s+([A-Z][a-z]?)\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})\s+(\-?[0-9]+\.[0-9]{3})'
#orb_regex = '\s+(([0-9]+)\-\s+)?([0-9]+)\s+[SPDFGH]'

spin_regex = '\s+(ALPHA[\+\-]BETA|BETA|ALPHA)\s+ELECTRONS'
matrix_regex =  '\s+([A-Z]+\s+MATRIX)\s+\-\s+CELL\s+N\.\s+([0-9]+)\(\s*(\-?[0-9]+)\s+(\-?[0-9]+)\s+(\-?[0-9]+)\s*\)'

NoSpin = 'nospin'
AplusB = 'a+b'
AminusB = 'a-b'



class CrystalMatrix:
	OVERLAP = 0
	DENSITY = 1
	NONE = -1

	def __init__(self):
#		self.cell = None
		self.matrix = None

#	@staticmethod
#	def read_lattice(lines):
#		for n in xrange(len(lines)):
#			if lattice_str in lines[n]:
#				break
#		vectors = []
#		for line in lines[n+1: n+4]:
#			ls = line.split()
#			assert len(ls) == 3
#			v = Vector([(float(x) * Units.BOHR / Units.UNIT) for x in ls])
#			vectors.append(v)
#		return vectors

#	@staticmethod
#	def read_cell_and_basis(lines, vm):
#		atoms = []
#		first_n = -1
#		last_n = -1
#		numorbmap = {}
#		for n in xrange(len(lines)):
#			if basis_str in lines[n]:
#				break
#		for line in lines[n+4:]:
#			m = re.match(atom_regex, line)
#			m1 = re.match(orb_regex, line)
#			if not line:
#				if atoms:
#					al = atoms[-1]
#					if last_n != -1 and first_n != -1:
#						numorb = last_n - first_n + 1
#						numorbmap[al.name()] = numorb
#					al.data()[AtomKeys.ORBITAL_COUNT] = numorbmap[al.name()]
#					al.data()[AtomKeys.FULL_VALENCE] = vm[al.name()]
#					al.data()[AtomKeys.ESTIMATED_VALENCE] = Shells.estimate_valence_byname(al.name())
#				break
#			elif m:
#				if atoms:
#					al = atoms[-1]
#					if last_n != -1 and first_n != -1:
#						numorb = last_n - first_n + 1
#						numorbmap[al.name()] = numorb
#					al.data()[AtomKeys.ORBITAL_COUNT] = numorbmap[al.name()]
#					al.data()[AtomKeys.FULL_VALENCE] = vm[al.name()]
#					al.data()[AtomKeys.ESTIMATED_VALENCE] = Shells.estimate_valence_byname(al.name())
#						
#				v = Vector([(float(x) * Units.BOHR / Units.UNIT) for x in [m.group(3), m.group(4), m.group(5)]])
#				name = m.group(2)
#				a = AtomVector(name, v)

#				atoms.append(a)
#				first_n = -1
#				last_n = -1
#			elif m1:
#				if first_n == -1:
#					if m1.group(1):
#						first_n = int(m1.group(2))
#					else:
#						first_n = int(m1.group(3))
#				last_n = int(m1.group(3))
#		return atoms

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
	def read_sdm_block(lines, n, atoms, satoms, prec):
		data = {}
		ls = [int(i) for i in lines[n].split()]
		atomsv = atoms
		atomsh = satoms
		h_indices = [CrystalMatrix.n_to_atomorb_tuple(i, atomsh) for i in ls]
		ntot = reduce((lambda x,y: x + y), [ca.data()[AtomKeys.ORBITAL_COUNT] for ca in atoms])
		v_indices = [CrystalMatrix.n_to_atomorb_tuple(i+1, atomsv) for i in range(ntot)]
		n += 1
		if lines[n] == '':
			n += 1
		for k in range(n, n + ntot):
			ls = lines[k].split()
			assert len(ls) == len(h_indices) + 1
			vi = k - n + 1
			assert int(ls[0]) == vi, vi
			for vj in range(1, len(ls)):
				value = float(ls[vj])
				if abs(value) > prec:
					avd = v_indices[vi - 1]
					ahd = h_indices[vj - 1]
					key = tuple(sorted([avd, ahd]))
					avd1 = ((avd[0][0], avd[0][1]-ahd[0][1], avd[0][2]-ahd[0][2], avd[0][3]-ahd[0][3]), avd[1])
					ahd1 = ((ahd[0][0], 0, 0, 0), ahd[1])
					key1 = tuple(sorted([avd1, ahd1]))
					data[key] = value
					data[key1] = value
		return data, k

	@staticmethod
	def read_sdm_cell(lines, n, cell, prec):
		m = re.match(matrix_regex, lines[n])
		assert m
		shifts = [int(m.group(i)) for i in (3,4,5)]
		
		if int(m.group(2)) % 25 == 0:
			logging.debug(lines[n])
		atoms = cell.cell
		satoms = cell.shift(shifts)
		n += 2
		datamap = {}
		while True:
			
			data, n = CrystalMatrix.read_sdm_block(lines, n, atoms, satoms, prec)
			datamap.update(data)
			n += 1
			if lines[n] != '':
				return datamap, n
			else:
				n += 1
				if re.match(matrix_regex, lines[n]) or re.match(spin_regex, lines[n]):
					return datamap, n

	@staticmethod
	def read_sdm(lines, n, cell, prec):
		datamap = {}
		while True:
			data, n = CrystalMatrix.read_sdm_cell(lines, n, cell, prec)
			datamap.update(data)
			if re.match(spin_regex, lines[n]):
				return datamap, n
			elif 'TTTTTTTTTTTTTT' in lines[n]:
				return datamap, n		

	@staticmethod
	def read_dm(lines, cell, prec):
		n = 0
		matrix_map = {}
		while True:
			mtype = NoSpin
			for k in range(n, len(lines)):
				m = re.match(spin_regex, lines[k])
				m1 = re.match(matrix_regex, lines[k])
				if m:
					logging.debug('****************')
					logging.debug(lines[k])
					logging.debug('****************')
					mt = m.group(1)
					if mt == 'ALPHA+BETA':
						mtype = AplusB
					elif mt == 'ALPHA-BETA':
						mtype = AminusB
					else:
						assert False, mt
					k1 = k + 1
					for k in range(k1, len(lines)):
						m1 = re.match(matrix_regex, lines[k])
						if m1:
							break
					break
				elif m1:
					mtype = NoSpin
					break
			data, n = CrystalMatrix.read_sdm(lines, k, cell, prec)
			assert mtype not in list(matrix_map.keys())
			matrix_map[mtype] = data
			if 'TTTTTTTTTTTTTT' in lines[n]:
				break
#		if matrix_map.keys() == NoSpin:
#			logging.debug(u'')
#			logging.debug(u'No spin matrix')
#			mat = {}
#			mat['a'] = {}
#			mat['a']['a'] = {}
#			mat['a']['a']['re'] = matrix_map[NoSpin]
#			mat
		if Counter(list(matrix_map.keys())) == Counter([AplusB, AminusB]):
			logging.debug('')
			logging.debug('2-spin matrix, converting')
			mat = {}
			mat['a'] = {}
			mat['b'] = {}
			mat['a']['a'] = {}
			mat['b']['b'] = {}
			mat['a']['a']['re'] = {}
			mat['b']['b']['re'] = {}
			keysp = set(matrix_map[AplusB].keys())
			keysm = set(matrix_map[AminusB].keys())
			for key in keysp&keysm:
				mat['a']['a']['re'][key] = (matrix_map[AplusB][key] + matrix_map[AminusB][key]) / 2
				mat['b']['b']['re'][key] = (matrix_map[AplusB][key] - matrix_map[AminusB][key]) / 2
			for key in keysp - keysm:
				mat['a']['a']['re'][key] = matrix_map[AplusB][key] / 2
				mat['b']['b']['re'][key] = matrix_map[AplusB][key] / 2
			for key in keysm - keysp:
				mat['a']['a']['re'][key] = matrix_map[AminusB][key] / 2
				mat['b']['b']['re'][key] = -matrix_map[AminusB][key] / 2
			
			return mat
		elif Counter(list(matrix_map.keys())) == Counter([AplusB]):
			logging.debug('')
			logging.debug('1-spin matrix, converting')
			mat = {}
			mat['a'] = {}
			mat['b'] = {}
			mat['a']['a'] = {}
			mat['b']['b'] = {}
			mat['a']['a']['re'] = {}
			mat['b']['b']['re'] = {}
			keysp = set(matrix_map[AplusB].keys())
			for key in keysp:
				mat['a']['a']['re'][key] = matrix_map[AplusB][key] / 2
				mat['b']['b']['re'][key] = matrix_map[AplusB][key] / 2
			return mat
		else:
			assert False, list(matrix_map.keys())

	@staticmethod
	def read_olp_block(lines, n, atoms, satoms, prec):
#		print lines[n]
		data = {}
		ls = [int(i) for i in lines[n].split()]
		atomsv = atoms
		atomsh = satoms
		h_indices = [CrystalMatrix.n_to_atomorb_tuple(i, atomsh) for i in ls]
		ntot = reduce((lambda x,y: x + y), [ca.data()[AtomKeys.ORBITAL_COUNT] for ca in atoms])
		nstart = ls[0]
		v_indices = [CrystalMatrix.n_to_atomorb_tuple(i+1, atomsv) for i in range(ntot)]
		n += 1
		if lines[n] == '':
			n += 1
		for k in range(n, n + ntot - nstart + 1):
			ls = lines[k].split()
			assert len(ls) == min(len(h_indices), k - n + 1) + 1
			vi = k - n + nstart
			assert int(ls[0]) == vi, vi
			for vj in range(1, len(ls)):
				value = float(ls[vj])
				if abs(value) > prec:
					avd = v_indices[vi - 1]
					ahd = h_indices[vj - 1]
					key = tuple(sorted([avd, ahd]))
					avd1 = ((avd[0][0], avd[0][1]-ahd[0][1], avd[0][2]-ahd[0][2], avd[0][3]-ahd[0][3]), avd[1])
					ahd1 = ((ahd[0][0], 0, 0, 0), ahd[1])
					key1 = tuple(sorted([avd1, ahd1]))
					data[key] = value
					data[key1] = value
		return data, k

	@staticmethod
	def read_olp_cell(lines, n, cell, prec):
		m = re.match(matrix_regex, lines[n])
		assert m
		shifts = [int(m.group(i)) for i in (3,4,5)]
		
		if int(m.group(2)) % 25 == 0:
			logging.debug(lines[n])
		atoms = cell.cell
		satoms = cell.shift(shifts)
		n += 2
		datamap = {}
		while True:
			data, n = CrystalMatrix.read_olp_block(lines, n, atoms, satoms, prec)
			datamap.update(data)
			n += 1
			if lines[n] != '':
				return datamap, n
			else:
				n += 1
				if re.match(matrix_regex, lines[n]):
					return datamap, n

	@staticmethod
	def read_olp(lines, cell, prec):
		n = 0
		matrix_map = {}
		for k in range(n, len(lines)):
			m1 = re.match(matrix_regex, lines[k])
			if m1:
				break
				datamap = {}
		n = k
		datamap = {}
		while True:
			data, n = CrystalMatrix.read_olp_cell(lines, n, cell, prec)
			datamap.update(data)
			if 'TTTTTTTTTTTTTT' in lines[n]:
				return datamap	

			
	@staticmethod
	def from_string(datastring, cell, mtype, prec):
		lines = datastring.splitlines()
#		vectors = CrystalMatrix.read_lattice(lines)
#		atoms = CrystalMatrix.read_cell_and_basis(lines, co.valence_map)
		cm = CrystalMatrix()
#		cm.cell = Cell(atoms, vectors)
		if mtype == CrystalMatrix.DENSITY:
			cm.matrix = CrystalMatrix.read_dm(lines, cell, prec)
		elif mtype == CrystalMatrix.OVERLAP:
			cm.matrix = CrystalMatrix.read_olp(lines, cell, prec)
		return cm

	@staticmethod
	def from_file(name, cell, mtype, prec=0):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Reading matrix from %s' % name)
		logging.info('*********************************************')
		logging.info('')
		with open(name) as f:
			return CrystalMatrix.from_string(f.read(), cell, mtype, prec)


