#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.util.xyz_format import read_xyz
from qcldm.util.units import Units
from qcldm.util.log_colorizer import init_log

IONIC_RADII = {
'O': 1.32,
'Y': 0.893,
'P': 0.35,
'Ca': 0.99,
'Nb': 0.69,
'U': 0.97,
'Y': 0.89,
}

#FORMULA_TEMPLATE = 'where(abs(x-%(R)f)>%(r)f,0,(1-((x-%(R)f)/%(r)f)**2)**0.5)'

FORMULA_TEMPLATE = '''
where(
  abs(x-%(R)f)>%(r)f,
  0,
  (%(r)f-absolute(x-%(R)f))/%(r)f*reciprocal(absolute(x-%(R)f)+1)/1
)
'''.replace("\n", "").replace(" ", "")

#FORMULA_TEMPLATE = '''
#where(
#  abs(x-%(R)f)>%(r)f,
#  0,
#  where(
#    x-%(R)f>0,
#    1-(1-((x-%(R)f-%(r)f)/%(r)f)**2)**0.5,
#    1-(1-((x-%(R)f+%(r)f)/%(r)f)**2)**0.5
#  )
#)
#'''.replace("\n", "").replace(" ", "")

init_log(sys.argv)

filename = sys.argv[1]
r0 = float(sys.argv[2]) * Units.BOHR

cluster = read_xyz(filename)

neighmap = {}
center = cluster.atoms[0]
for atom in cluster.atoms[0:]:
	r = atom.distance(center)
	if r0 < r:
		continue
	if atom.name() not in neighmap.keys():
		neighmap[atom.name()] = {}
	foundr = False
	for rr in neighmap[atom.name()].keys():
		if abs(rr - r) < 0.01:
			foundr = True
			neighmap[atom.name()][rr] += ' + ' + FORMULA_TEMPLATE % {'R' : atom.distance(center) / Units.BOHR, 'r' : IONIC_RADII[atom.name()] / Units.BOHR}
			break
	if not foundr:
		neighmap[atom.name()][r] = FORMULA_TEMPLATE % {'R' : atom.distance(center) / Units.BOHR, 'r' : IONIC_RADII[atom.name()] / Units.BOHR}
			

with open('cluster_neghbour_formulas', 'w') as o:
	for k in sorted(neighmap.keys()):
		o.write("%s\n" % (k))
		for r in sorted(neighmap[k].keys()):
			o.write("%f\n%s\n" % (r / Units.BOHR, neighmap[k][r]))
