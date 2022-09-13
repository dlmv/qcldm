#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.cube_format.gaussian_cube import GaussianCube
from qcldm.cube_format.cube_operations import transformed_copy, rescale_simple, subtract
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from math3d import Vector

from qcldm.crystal_format.crystal_out import CrystalOut
from qcldm.crystal_format.crystal_meta import CrystalMeta

init_log(sys.argv)

c = CrystalMeta()
c.load('.')

co = CrystalOut.from_file(c.out_file)
write_xyz(co.cell.cell, 'cell.xyz')
write_xyz(co.cell.supercell, 'supercell.xyz')

cube = GaussianCube.from_file('crystal.cube')

cube.celltype = GaussianCube.TYPE_PERIODIC_WITH_BORDER

a1 = Vector([0.00000000,1.29001007,-1.29858367])
a2 = Vector([0.00000000,-1.29001007,-3.89575101])

a1 = co.cell.validate_point(a1)
a2 = co.cell.validate_point(a2)

symop_to_use = None

for symop in co.cell.symops:
	na1 = co.cell.apply_symop_to_point(a1, symop, True)
	if na1 == a2:
		symop_to_use = symop
		print(symop)
		break

trans = transformed_copy(cube, co.cell, symop_to_use)

trans.to_file('transformed.cube')

cube.celltype = GaussianCube.TYPE_CLUSTER

resc = rescale_simple(trans, cube)

trans.to_file('rescaled.cube')

diff = subtract(resc, cube)


diff.to_file('diff.cube')

