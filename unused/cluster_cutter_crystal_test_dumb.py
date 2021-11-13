#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.crystal_format.crystal_matrix import CrystalMatrix
from qcldm.crystal_format.crystal_out import CrystalOut
from qcldm.crystal_format.crystal_meta import CrystalMeta
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.structures.cluster_embedding import Cluster
from qcldm.structures.bader_reader import read_baders
from qcldm.structures.atom_vector import  AtomKeys

init_log(sys.argv)

c = CrystalMeta()
c.load('.')

co = CrystalOut.from_file(c.out_file)
#ocm = CrystalMatrix.from_file('overlap_ypo4.outp', co.cell, CrystalMatrix.NONE)
write_xyz(co.cell.cell, 'cell.xyz')
#write_xyz(co.cell.supercell, 'supercell.xyz')

read_baders(co.cell)

num = int(sys.argv[1])
layers = int(sys.argv[2])
electro = int(sys.argv[3])
center = co.cell.cell[num - 1]
centers = [center]

tmpshells = co.cell.neighbours.neighbours_cluster(centers, 2)
for a in tmpshells[-1]:
	if a.name() != center.name():
		centers.append(a)

cluster = Cluster(co.cell, centers, layers, electro)

dirname = "cluster_test_dumb%d_%d_%d" % (num, layers, electro)

key = AtomKeys.BADER_CHARGE

cluster.estimate_charges_dumb(key)

cluster.write_structure(dirname)
cluster.write_charges(key, dirname)
cluster.write_embedding(key, dirname)



