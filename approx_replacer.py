#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.util.log_colorizer import init_log
from qcldm.openmx_format.pao_format import PAO
from qcldm.gauss_format.gauss_format import GaussFormat

init_log(sys.argv)

paof = sys.argv[1]
gauf = paof.replace(".pao", ".gauss")
appf = paof.replace(".pao", "A.pao")


basis = PAO.from_file(paof)
gs = GaussFormat.from_gaussian94(open(gauf).readlines())

fs = []
ns = {}

for g in gs:
	l = g.fs[0][1].l
	if l not in ns.keys():
		ns[l] = 1
	f = g.to_numeric(basis.grid)
	f.l = l
	f.n = ns[l]
	ns[l] += 1
	fs.append(f)

basis.set_functions(fs, basis.grid)

PAO.to_file(basis, appf)

















