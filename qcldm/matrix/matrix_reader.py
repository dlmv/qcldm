import re, os, sys, logging

from ..structures.cell import CellAtom

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

def read_cell_atom(line, c):
	ls = filter(None, re.split('\s*', line))
	a = CellAtom(c, int(ls[1]), [int(ls[2]), int(ls[3]), int(ls[4])])
	return a

def read_matrix_atoms(name, c):
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Reading matrix atoms')
	logging.info(u'*********************************************')
	logging.info(u'')
	with open(name, 'r') as f:
		atoms = []
		lines = f.readlines()
		n = int(lines[0][len(HEADER):])
		for line in lines[1:n+1]:
			a = read_cell_atom(line, c)
			atoms.append(a)
		logging.debug(u'Total atoms: %d' % len(atoms))
		return atoms

def read_header(line, atoms):
	orbs = []
	ls = filter(None, re.split('[\s*|\**]', line))
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
		ls = filter(None, re.split('\s*', lines[n]))
		atom1 = CellAtom(c, int(ls[0]), [0, 0, 0])
		o1 = (atom1.tuple_data(), int(ls[1]))
		for i in xrange(len(ls[2:])):
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
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Reading matrix: %s' % name)
	logging.info(u'*********************************************')
	logging.info(u'')
	with open(name, 'r') as f:
		datamap = {}
		lines = f.readlines()
		read_header(lines[4], mat_atoms)
		n = 1
		while n < len(lines):
			data, n = read_matrix_part(lines, n, c, mat_atoms, prec)
			logging.debug(u'Progress: %d/%d' % (n, len(lines)))
			datamap.update(data)

	return datamap




