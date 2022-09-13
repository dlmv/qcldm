#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)


c = GaussianCube.from_file(sys.argv[1])
#c.rescale(Vector(0,0,0), [s/2 for s in c.size], c.vectors)
nvs = [v.copy() for v in c.vectors]
for nv, v in zip(nvs, c.vectors):
	for v1 in c.vectors:
		if v1 != v:
			nv += v1/2
print(nvs)
c.rescale(c.origin, c.size, nvs)
#c.rescale(c.origin, c.size, c.vectors)
c.to_file('cube.cube')



