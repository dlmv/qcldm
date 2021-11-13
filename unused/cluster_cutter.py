#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.openmx_format.dat_format import DAT_INPUT
from qcldm.applications.openmx_rewriter import rewrite_files
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.matrix.matrix_reader import read_matrices
from qcldm.structures.cluster_embedding import Cluster
from qcldm.structures.bader_reader import read_baders

init_log(sys.argv)

d = DAT_INPUT.from_file('temporal_12345.input')
write_xyz(d.cell.cell, 'cell.xyz')

read_baders(d.cell)

dm, olp, atoms = read_matrices(d.cell, 1e-3)

num = int(sys.argv[1])
layers = int(sys.argv[2])
electro = int(sys.argv[3])

centers = [d.cell.cell[num - 1]]

cluster = Cluster(d.cell, centers, layers, electro)
cluster.estimate_charges(dm, olp)

dirname = "cluster%d_%d_%d" % (num, layers, electro)

cluster.write_structure(dirname)
cluster.write_charges(dirname)
cluster.write_embedding(dirname)

rewrite_files(d, cluster, dirname)




