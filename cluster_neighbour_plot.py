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

#FORMULA_TEMPLATE = 'where(abs(x-%(R)f)>%(r)f,0,(1-((x-%(R)f)/%(r)f)**2)**0.5)'

FORMULA_TEMPLATE = '''
where(
  abs(x-%(R)f)>%(r)f,
  0,
  where(
    x-%(R)f>0,
    1-(1-((x-%(R)f-%(r)f)/%(r)f)**2)**0.5,
    1-(1-((x-%(R)f+%(r)f)/%(r)f)**2)**0.5
  )
)
'''.replace("\n", "").replace(" ", "")

init_log(sys.argv)

filename = sys.argv[1]
r = float(sys.argv[2]) * Units.BOHR

cluster = read_xyz(filename)

neighmap = {}
center = cluster.atoms[0]
for atom in cluster.atoms[0:]:
	if r < atom.distance(center):
		continue
	if atom.name() not in neighmap.keys():
		neighmap[atom.name()] = ''
	else:
		neighmap[atom.name()] += ' + '
	neighmap[atom.name()] += FORMULA_TEMPLATE % {'R' : atom.distance(center) / Units.BOHR, 'r' : IONIC_RADII[atom.name()] / Units.BOHR / 2}

with open('cluster_neghbour_formulas', 'w') as o:
	for k in sorted(neighmap.keys()):
		o.write("%s\n%s\n" % (k, neighmap[k]))
