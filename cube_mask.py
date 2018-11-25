#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import masked_atomsphere
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)

r = float(sys.argv[1])

cube = GaussianCube.from_file('diff.cube')

c = masked_atomsphere(cube, 1, r)


c.to_file('masked.cube')



