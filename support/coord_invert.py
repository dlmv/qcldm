#!/usr/bin/python 
import os, sys, re



filename = sys.argv[1]
icoords = [float(x) / 0.529177249 for x in sys.argv[2:5]]

with open(filename) as coord:
	with open(filename + '_inverted', 'w') as coord_inv:
		lines = coord.read().splitlines()
		assert lines[0] == '$coord'
		assert lines[-1] == '$end'
		coord_inv.write('$coord\n')
		for line in lines[1:-1]:
			ls = line.split()
			x, y, z = [i*2 - float(r) for r, i in zip(ls[:3], icoords)]
			name = ls[3]
			line1 = '%16.10f %16.10f% 16.10f %2s \n' % (x, y, z, name)
			coord_inv.write(line1)
		
		coord_inv.write('$end\n')























