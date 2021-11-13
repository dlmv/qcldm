#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import subtract
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)

cube1 = GaussianCube.from_file('rescaled.cube')
cube2 = GaussianCube.from_file('cluster.cube')

#cs = GaussianCube.from_file('mag.cube')

c = subtract(cube1, cube2)


c.to_file('diff.cube')



