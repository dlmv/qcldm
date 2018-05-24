#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.crystal_format.crystal_matrix import CrystalMatrix
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.structures.cluster_embedding import Cluster
from qcldm.structures.bader_reader import read_baders

init_log(sys.argv)

cm = CrystalMatrix.from_file('density_ypo4.outp', 1e-2)
write_xyz(cm.cell.cell, 'cell.xyz')
write_xyz(cm.cell.supercell, 'supercell.xyz')

#read_baders(d.cell)

#dm, olp, atoms = read_matrices(d.cell, 1e-3)

#num = int(sys.argv[1])
#layers = int(sys.argv[2])
#electro = int(sys.argv[3])

#centers = [d.cell.cell[num - 1]]

#cluster = Cluster(d.cell, centers, layers, electro)
#cluster.estimate_charges(dm, olp)

#dirname = "cluster%d_%d_%d" % (num, layers, electro)

#cluster.write_structure(dirname)
#cluster.write_charges(dirname)
#cluster.write_embedding(dirname)

#rewrite_files(d, cluster, dirname)




