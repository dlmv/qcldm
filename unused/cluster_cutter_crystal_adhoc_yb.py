#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.crystal_format.crystal_out import CrystalOut
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.structures.cluster_embedding import Cluster
from qcldm.structures.atom_vector import  AtomKeys

init_log(sys.argv)

co = CrystalOut.from_file('out')
for a in co.cell.atoms:
	if a.name() == 'Yb':
		a.data()[AtomKeys.ESTIMATED_VALENCE] = 3
		a.data()[AtomKeys.ESTIMATED_CHARGE] = 3
write_xyz(co.cell.cell, 'cell.xyz')
write_xyz(co.cell.supercell, 'supercell.xyz')

num = int(sys.argv[1])
layers = int(sys.argv[2])
electro = int(sys.argv[3])

centers = [co.cell.cell[num - 1]]

cluster = Cluster(co.cell, centers, layers, electro)

key = AtomKeys.MULLIKEN_CHARGE

cluster.estimate_charges_dumb(key)

dirname = "cluster_%d_%d_%d" % (num, layers, electro)

cluster.write_structure(dirname)
cluster.write_charges(key, dirname)
cluster.write_embedding(key, dirname)





