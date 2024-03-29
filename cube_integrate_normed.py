#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import integrate_in_sphere_range_normed
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)

filename = sys.argv[1]
r = float(sys.argv[3])
atomnum = int(sys.argv[2])

cube = GaussianCube.from_file(filename)

rs, res = integrate_in_sphere_range_normed(cube, atomnum, r, 0.1)
with open(filename + ".integrated", "w") as o:
	for r, v in zip(rs, res):
		o.write("%8.2f %12.6e\n" % (r, v))






