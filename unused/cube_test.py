#!/usr/bin/python
import re, sys, math, logging
import numpy as np
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import rescale_simple, rescale_medium
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

init_log(sys.argv)


source = GaussianCube.from_file('crystal.cube')

source.celltype = GaussianCube.TYPE_PERIODIC_WITH_BORDER

#print source.size

#for x in range(source.size[0]):
#	for y in range(source.size[1]):
#		for z in range(source.size[2]):
#			if source.data[x,y,z] > 100:
#				print x, y, z, source.voxel_cuboid([x,y,z]).center

p1 = [0.000000, 2.437766, 7.361902]

p2 = [0.000000, 2.437766, -2.453967]

print((source.voxel_cuboid(source.find_voxel(p1)).center - np.array(p1)) - (source.voxel_cuboid(source.find_voxel(p2)).center - np.array(p2)))

#print source.find_voxel(p1)
#print source.voxel_cuboid(source.find_voxel(p1)).center
#print source.point_value(p1)

#print source.point_value()


