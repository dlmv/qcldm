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
write_xyz(co.cell.cell, 'cell.xyz')
#write_xyz(cm.cell.supercell, 'supercell.xyz')

print(co.cell.cell)

num = 1

centers = [co.cell.cell[num - 1]]

tmpshells = co.cell.neighbours.neighbours_cluster(centers, 1)

atoms = []

for t in tmpshells:
	atoms.extend(t)

write_xyz(atoms, 'tt.xyz')





