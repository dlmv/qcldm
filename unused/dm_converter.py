#!/usr/bin/python
import re, sys, os, logging
sys.dont_write_bytecode = True

from qcldm.openmx_format.dat_format import DAT_INPUT
from qcldm.openmx_format.openmx_orbital_order import Openmx_Order
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.matrix.matrix_reader import read_matrices, merge_matrices
from qcldm.applications.complex_matrix import convert_atom_matrix

init_log(sys.argv)

d = DAT_INPUT.from_file('temporal_12345.input')
write_xyz(d.cell.cell, 'cell.xyz')

if len(d.cell.cell) != 1:
	logging.error(u'Only for a single atom!!!')
	sys.exit(1)

dms, olp, atoms = read_matrices(d.cell)

DM, OLP = merge_matrices(dms, olp, atoms)


convert_atom_matrix(DM, OLP, d.cell.cell[0], Openmx_Order())






