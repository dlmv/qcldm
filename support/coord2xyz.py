#!/usr/bin/python 
import os, sys, re



filename = sys.argv[1]

with open(filename) as coord:
	with open(filename + '.xyz', 'w') as xyz:
		lines = coord.read().splitlines()
		assert lines[0] == '$coord'
		assert lines[-1] == '$end'
		for line in lines[1:-1]:
			ls = line.split()
			x, y, z = [float(r) * 0.529177249 for r in ls[:3]]
			name = ls[3]
			if name == 'zz' or name == 'ZZ' or name == 'q' or name == 'Q':
				name = 'X'
			line1 = '%2s %16.10f %16.10f% 16.10f\n' % (name, x, y, z)
			xyz.write(line1)























