#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz

init_log(sys.argv)


c = GaussianCube.from_file(sys.argv[1])
c.rescale(c.origin, c.size, c.vectors)
c.to_file('cube.cube')



