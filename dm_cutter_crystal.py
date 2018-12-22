#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.crystal_format.crystal_matrix import CrystalMatrix
from qcldm.crystal_format.crystal_out import CrystalOut
from qcldm.crystal_format.crystal_meta import CrystalMeta
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.applications.oneprop import prepare_oneprop_crystal

init_log(sys.argv)

num = int(sys.argv[1])
rc = float(sys.argv[2])
prec = float(sys.argv[3])

c = CrystalMeta()
c.load('.')

co = CrystalOut.from_file(c.out_file)
co.get_cutoffs(prec)

dcm = CrystalMatrix.from_file(c.dm_file, co.cell, CrystalMatrix.DENSITY, 0)
write_xyz(co.cell.cell, 'cell.xyz')


centers = [co.cell.cell[num - 1]]
cluster = co.cell.neighbours.neighbours_cluster(centers, 6, {})
atoms = []
for shell in cluster:
	atoms.extend(shell)
prepare_oneprop_crystal(co, atoms, dcm.matrix, num, rc)







