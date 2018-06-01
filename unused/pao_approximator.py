#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.util.log_colorizer import init_log
from qcldm.atom.shells import Shells
from qcldm.openmx_format.pao_format import PAO
from qcldm.applications.gauss_approximator import approximate

init_log(sys.argv)

paof = sys.argv[1]
gauf = paof.replace(".pao", ".gauss")
rgauf = paof.replace(".pao", ".gaussr")
p = PAO.from_file(paof)

with open(gauf, 'w') as fl:
	fl.write("X  0\n")
with open(rgauf, 'w') as fl:
	fl.write("X  0\n")

writefile = '-f' in sys.argv

nmin, nmax = [int(n) for n in sys.argv[2].split('-')]

for l in sorted(p.functions.keys()):
	for n in sorted(p.functions[l].keys()):
		f = p.functions[l][n]
		cg, res = approximate(f, 1e-3, 16, 50000, nmin, nmax, writefile)
		with open(gauf, 'a') as fl:
			fl.write("%s  %2d 1.00   \n" % (Shells.SHELLS[l], len(cg.fs)))
			for n, gf in cg.fs:
				fl.write("    %25.18f    %25.18f\n" % (gf.a, n))

		with open(rgauf, 'a') as fl:
			fl.write("%s  %2d 1.00                                                  R = %e, NORM = %f\n" % (Shells.SHELLS[l], len(cg.fs), res, cg.norm()))
			for n, gf in cg.fs:
				fl.write("    %25.18f    %25.18f\n" % (gf.a, n))















