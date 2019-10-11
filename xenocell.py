#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.crystal_format.crystal_out import CrystalOut
from qcldm.crystal_format.crystal_meta import CrystalMeta
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.embedding.cluster import Cluster
from qcldm.embedding.embedding_settings import EmbeddingSettings
from qcldm.structures.bader_reader import read_baders
from qcldm.structures.atom_vector import AtomKeys

init_log(sys.argv)

c = CrystalMeta()
c.load('.')

co = CrystalOut.from_file(c.out_file)
write_xyz(co.cell.cell, 'cell.xyz')
write_xyz(co.cell.supercell, 'supercell.xyz')

atoms = []
vectors = [6.91898051, 6.91898051, 6.06871692]
zero = [-2.3, -3, -2.3]
for a in co.cell.extended_cell(3):
	ok = True
	for k in range(3):
		if a.position()[k] < zero[k] or a.position()[k] >= zero[k] + vectors[k]:
			ok = False
			break
	if ok:
		atoms.append(a)
write_xyz(atoms, 'ortho1.xyz')
	




