#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz

init_log(sys.argv)


c = GaussianCube.from_file('cube.cube')
#c.scale([3,3,3])
c.multiply([3,3,3])
c.to_file('cube1.cube')



