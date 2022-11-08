import os, sys, re, shutil


cube_template = '''$%sdensity-cube
x0= %.6f y0= %.6f z0= %.6f
nx= %d ny= %d nz= %d
dx= %.6f dy= %.6f dz= %.6f
'''

def create_cube_input(coords, r, spin):
	start = coords
	vectors = [0.05] * 3
	x0, y0, z0 = [int((coords[i] - r - start[i]) / vectors[i]) * vectors[i] + start[i] for i in range(3)]
	x1, y1, z1 = [int((coords[i] + r - start[i]) / vectors[i]) * vectors[i] + start[i] for i in range(3)]

	dx, dy, dz = vectors
	nx, ny, nz = [int(([x1, y1, z1][i]-[x0, y0, z0][i])/vectors[i]) for i in range(3)]
	res = cube_template % ('spin' if spin else '', x0, y0, z0, nx, ny, nz, dx, dy, dz)
	return res

def read_center_coords():
	with open('coord') as crf:
		head = [next(crf) for x in range(2)]
		coords = [float(x) for x in head[1].split()[:-1]]
		return coords

def rewrite_control(cube_input):
	skipping = False
	written = False
	shutil.copy('control', 'control.bak')
	with open('control.bak') as inp:
		with open('control', 'w') as outp:
			for line in inp:
				if skipping and line[0] == '$':
					skipping = False
				if line.startswith('$dft-functional') and not written:
					line += cube_input
					written = True
				elif line.startswith('$spindensity-cube') or  line.startswith('$density-cube'):
					skipping = True
				if not skipping:
					outp.write(line)
			
		
r = float(sys.argv[1])
spin = len(sys.argv) > 2 and sys.argv[2] == 'spin'
			
coords =read_center_coords()

cube_input = create_cube_input(coords, r, spin)
#print(cube_input)

rewrite_control(cube_input)


















