#!/usr/bin/python
import re, sys, os
sys.dont_write_bytecode = True

from qcldm.openmx_format.dat_format import DAT_INPUT
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.applications.oneprop import prepare_oneprop
from qcldm.matrix.matrix_reader import read_matrices

init_log(sys.argv)

d = DAT_INPUT.from_file('temporal_12345.input')
write_xyz(d.cell.cell, 'cell.xyz')

dm, olp, atoms = read_matrices(d.cell)

num = int(sys.argv[1])
rc = float(sys.argv[2])

prepare_oneprop(d, atoms, dm, num, rc)





