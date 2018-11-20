#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import rescale_simple
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)


target = GaussianCube.from_file('dens.cube')

source = GaussianCube.from_file('cube_fersmite_dat.DENS_CUBE')

source.celltype = GaussianCube.TYPE_PERIODIC_WITH_BORDER

#cs = GaussianCube.from_file('mag.cube')

c = rescale_simple(source, target)


c.to_file('rescaled.cube')



