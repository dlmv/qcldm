#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.openmx_format.dat_format import DAT_INPUT
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.structures.cluster_embedding import Cluster

init_log(sys.argv)

d = DAT_INPUT.from_file('test.dat')
write_xyz(d.cell.cell, 'cell.xyz')

num = int(sys.argv[1])
layers = int(sys.argv[2])
electro = int(sys.argv[3])
center = d.cell.cell[num - 1]
centers = [center]

tmpshells = d.cell.neighbours.neighbours_cluster(centers, 2)
for a in tmpshells[-1]:
	if a.name() != center.name():
		centers.append(a)

cluster = Cluster(d.cell, centers, layers, electro)

dirname = "cluster_test%d_%d_%d" % (num, layers, electro)

cluster.load_atoms_dummy()
cluster.write_structure(dirname)
