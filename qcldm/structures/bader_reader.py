import re
from ..structures.atom_vector import AtomKeys

def read_baders(cell, name='ACF.dat'):
	with open(name, 'r') as f:
		lines = f.readlines()
		bader_occs = []
		n = 0
		while '------------------' not in lines[n]:
			n += 1
		n += 1
		while '------------------' not in lines[n]:
			ls = filter(None, re.split('\s*', lines[n]))
			bc =float(ls[4])
			bader_occs.append(bc)
			n += 1
	for b, a in zip(bader_occs, cell.atoms):
		bc = a.data()[AtomKeys.FULL_VALENCE] - b
		a.data()[AtomKeys.BADER_CHARGE] = bc
