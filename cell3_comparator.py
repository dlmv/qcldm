#!/usr/bin/python
import sys
sys.dont_write_bytecode = True

#from qcldm.openmx_format.dat_format import DAT_INPUT
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz, read_xyz
from qcldm.structures.cell_difference import compare_cells3

init_log(sys.argv)

name1 = sys.argv[1]
name2 = sys.argv[2]
name3 = sys.argv[3]


cell1 = read_xyz(name1)
cell2 = read_xyz(name2)
cell3 = read_xyz(name3)

compare_cells3(cell1, cell2, cell3)





