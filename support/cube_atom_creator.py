import os, sys, re, shutil

cube_template = '''$density-cube
x0= %.6f y0= %.6f z0= %.6f
nx= %d ny= %d nz= %d
dx= %.6f dy= %.6f dz= %.6f
'''

def create_cube_input(vectors, coords, r):
	x0, y0, z0 = [int(-r / vectors[i]) * vectors[i] + coords[i] for i in range(3)]
	x1, y1, z1 = [int(r / vectors[i]) * vectors[i] + coords[i] for i in range(3)]
	dx, dy, dz = vectors
	nx, ny, nz = [int(([x1, y1, z1][i]-[x0, y0, z0][i])/vectors[i]) for i in range(3)]
	res = cube_template % (x0, y0, z0, nx, ny, nz, dx, dy, dz)
	return res

def read_center_coords(i):
	with open('coord') as crf:
		line = crf.read().splitlines()[i]
		coords = [float(x) for x in line.split()[:-1]]
		return coords

def rewrite_control(cube_input):
	skipping = False
	written = False
	shutil.copy('control', 'control.bak')
	with open('control.bak') as inp:
		with open('control', 'w') as outp:
			for line in inp.xreadlines():
				if skipping and line[0] == '$':
					skipping = False
				if line.startswith('$dft-functional') and not written:
					line += cube_input
					written = True
				elif line.startswith('$density-cube'):
					skipping = True
				if not skipping:
					outp.write(line)
			
		
i, r = int(sys.argv[1]), float(sys.argv[2])

vectors = [0.05,0.05,0.05]
			
coords = read_center_coords(i)
print coords

cube_input = create_cube_input(vectors, coords, r)

rewrite_control(cube_input)


















