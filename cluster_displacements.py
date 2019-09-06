#!/usr/bin/python
import sys
sys.dont_write_bytecode = True

from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz, read_xyz
from qcldm.structures.cell_difference import get_dispacements

init_log(sys.argv)

name1 = sys.argv[1]
name2 = sys.argv[2]

cell1 = read_xyz(name1)
cell2 = read_xyz(name2)

sd, mxd, mxa = get_dispacements(cell1, cell2)

print "Standard deviation: %f" % sd
print "Max deviation: %f" % mxd
print mxa




