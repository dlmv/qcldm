#!/usr/bin/python
import sys
sys.dont_write_bytecode = True

#from qcldm.openmx_format.dat_format import DAT_INPUT
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz, read_xyz
from qcldm.structures.cell_difference import compare_cells

init_log(sys.argv)

name1 = sys.argv[1]
name2 = sys.argv[2]

#dat1 = DAT_INPUT.read_from_file(name1)
#dat2 = DAT_INPUT.read_from_file(name2)

cell1 = read_xyz(name1)
cell2 = read_xyz(name2)

sd, mxd, a, nb = compare_cells(cell1, cell2)

print("Standard deviation: %f" % sd)
print("Max deviation: %f" % mxd)
print(a)
print(nb)




