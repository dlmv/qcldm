#!/usr/bin/python
import re, sys, os, logging
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.turbomole_format.turbo_orbital_order import Turbo_Order
from qcldm.applications.complex_matrix import convert_atom_matrix
from qcldm.atom.shells import Shells
from qcldm.util.log_colorizer import init_log

init_log(sys.argv)

c = ControlFormat.from_file('control')

if len(c.cell.cell) != 1:
	logging.error(u'Only for a single atom!!!')
	sys.exit(1)
b = c.bases.values()[0]
olp = b.internal_overlap()
oar = b.orbarray()
OLP = []
i1 = 0
for l1, nl1 in enumerate(b.orbarray()):
	for n1 in range(nl1 / (2 * l1 + 1)):
		for m1 in range(2 * l1 + 1):
			for s1 in range(2):
				olpline = []
				i2 = 0
				for l2, nl2 in enumerate(b.orbarray()):
					for n2 in range(nl2 / (2 * l2 + 1)):
						for m2 in range(2 * l2 + 1):
							for s2 in range(2):
								se = 0
								if m1 == m2 and l1 == l2 and s1 == s2:
									se = b.internal_overlap()[i1][i2]
								olpline.append(se)
						i2 += 1
				OLP.append(olpline)
		i1 += 1

convert_atom_matrix(c.dm, OLP, c.cell.cell[0], Turbo_Order())


