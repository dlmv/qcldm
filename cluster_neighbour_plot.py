#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.util.xyz_format import read_xyz
from qcldm.util.units import Units
from qcldm.util.log_colorizer import init_log

IONIC_RADII = {
'O': 1.32,
'Y': 0.893,
'P': 0.35
}

FORMULA_TEMPLATE = 'exp(-((x-%f)**2)/(2*%f**2))'

init_log(sys.argv)

filename = sys.argv[1]
r = float(sys.argv[2]) * Units.BOHR

cluster = read_xyz(filename)

neighmap = {}
center = cluster.atoms[0]
for atom in cluster.atoms[1:]:
	if r < atom.distance(center):
		continue
	if atom.name() not in neighmap.keys():
		neighmap[atom.name()] = ''
	else:
		neighmap[atom.name()] += ' + '
	neighmap[atom.name()] += FORMULA_TEMPLATE % (atom.distance(center) / Units.BOHR, IONIC_RADII[atom.name()] / Units.BOHR / 10)

with open('cluster_neghbour_formulas', 'w') as o:
	for k in sorted(neighmap.keys()):
		o.write("%s\n%s\n" % (k, neighmap[k]))
