#!/usr/bin/python
import re, sys, os, logging
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.turbomole_format.turbo_orbital_order import Turbo_orbital_order
from qcldm.applications.complex_matrix import convert_atom_matrix
from qcldm.util.log_colorizer import init_log

init_log(sys.argv)

d = ControlFormat.from_file('control')

if len(d.cell.cell) != 1:
	logging.error(u'Only for a single atom!!!')
	sys.exit(1)


convert_atom_matrix(d.dm, d.cell.cell[0], Turbo_orbital_order())


