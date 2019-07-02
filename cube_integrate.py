#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import integrate_in_sphere
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)

filename = sys.argv[1]
r = float(sys.argv[2])

cube = GaussianCube.from_file(filename)

n = integrate_in_sphere(cube, 1, r)
print n




