import re, os, sys, logging

from ..structures.cell import CellAtom
from ..structures.atom_vector import AtomKeys

HEADER = 'all atoms: '


def read_matrices(cell, precision=0):
	matrix_atoms = read_matrix_atoms('atoms.mat', cell)
	dm = {}
	dm['a'] = {}
	dm['b'] = {}
	dm['a']['a'] = {}
	dm['a']['b'] = {}
	dm['b']['b'] = {}
	olp = read_matrix('olp.mat', cell, matrix_atoms, precision)
	dm['a']['a']['re'] = read_matrix('dm_aa_re.mat', cell, matrix_atoms, precision)
	dm['b']['b']['re'] = read_matrix('dm_bb_re.mat', cell, matrix_atoms, precision) or dm['a']['a']['re']
	dm['a']['b']['re'] = read_matrix('dm_ab_re.mat', cell, matrix_atoms, precision) or {}
	dm['a']['b']['im'] = read_matrix('dm_ab_im.mat', cell, matrix_atoms, precision) or {}
	dm['a']['a']['im'] = read_matrix('dm_aa_im.mat', cell, matrix_atoms, precision) or {}
	dm['b']['b']['im'] = read_matrix('dm_bb_im.mat', cell, matrix_atoms, precision) or {}
	return dm, olp, matrix_atoms

def merge_matrices(dms, olp, matrix_atoms):
	on = 0
	for a in matrix_atoms:
		on += a.data()[AtomKeys.ORBITAL_COUNT]
	DM = [[0] * on*2 for x in range(on*2)]
	OLP = [[0] * on*2 for x in range(on*2)]
	i = 0
	for a1 in matrix_atoms:
		for o1 in range(a1.data()[AtomKeys.ORBITAL_COUNT]):
			j = 0
			for a2 in matrix_atoms:
				for o2 in range(a2.data()[AtomKeys.ORBITAL_COUNT]):
					key = matrix_key(a1, a2, o1, o2)
					DM[2*i][2*j] = (dms['a']['a']['re'].get(key) or 0) + 1j * (dms['a']['a']['im'].get(key) or 0)
					DM[2*i + 1][2*j] = (dms['a']['b']['re'].get(key) or 0) + 1j * (dms['a']['b']['im'].get(key) or 0)
					DM[2*i][2*j + 1] = (dms['a']['b']['re'].get(key) or 0) - 1j * (dms['a']['b']['im'].get(key) or 0)
					DM[2*i + 1][2*j + 1] = (dms['b']['b']['re'].get(key) or 0) + 1j * (dms['b']['b']['im'].get(key) or 0)
					OLP[2*i][2*j] = (olp.get(key) or 0)
					OLP[2*i + 1][2*j + 1] = (olp.get(key) or 0)
					j += 1
			i += 1
	return DM, OLP

def matrix_key(a1, a2, o1, o2):
	return tuple(sorted([(a1.tuple_data(), o1 + 1), (a2.tuple_data(), o2 + 1)]))

def read_cell_atom(line, c):
	ls = [_f for _f in re.split('\s*', line) if _f]
	a = CellAtom(c, int(ls[1]), [int(ls[2]), int(ls[3]), int(ls[4])])
	return a

def read_matrix_atoms(name, c):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Reading matrix atoms')
	logging.info('*********************************************')
	logging.info('')
	with open(name, 'r') as f:
		atoms = []
		lines = f.readlines()
		n = int(lines[0][len(HEADER):])
		for line in lines[1:n+1]:
			a = read_cell_atom(line, c)
			atoms.append(a)
		logging.debug('Total atoms: %d' % len(atoms))
		return atoms

def read_header(line, atoms):
	orbs = []
	ls = [_f for _f in re.split('[\s*|\**]', line) if _f]
	for an, on in [(int(ls[i]), int(ls[i + 1])) for i in range(0, len(ls), 2)]:
		atom = atoms[an - 1]
		o = (atom.tuple_data(), on)
		orbs.append(o)
	return orbs

def read_matrix_part(lines, n, c, atoms, prec):
	while("****" not in lines[n]):
		n += 1
	orbs = read_header(lines[n], atoms)
	n += 1
	data = {}
	while n < len(lines) and 'matrix' not in lines[n]:
		ls = [_f for _f in re.split('\s*', lines[n]) if _f]
		atom1 = CellAtom(c, int(ls[0]), [0, 0, 0])
		o1 = (atom1.tuple_data(), int(ls[1]))
		for i in range(len(ls[2:])):
			o2 = orbs[i]
			key = tuple(sorted([o1, o2]))
			value = float(ls[i + 2])
			if abs(value) > prec:
				data[key] = value
		n += 1
	return data, n

def read_matrix(name, c, mat_atoms, prec=0):
	if not os.path.exists(name) or not mat_atoms:
		return None
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Reading matrix: %s' % name)
	logging.info('*********************************************')
	logging.info('')
	with open(name, 'r') as f:
		datamap = {}
		lines = f.readlines()
		read_header(lines[4], mat_atoms)
		n = 1
		while n < len(lines):
			data, n = read_matrix_part(lines, n, c, mat_atoms, prec)
			logging.debug('Progress: %d/%d' % (n, len(lines)))
			datamap.update(data)

	return datamap




