import os, sys, re, shutil

cube_template = '''$density-cube
x0= %.6f y0= %.6f z0= %.6f
nx= %d ny= %d nz= %d
dx= %.6f dy= %.6f dz= %.6f
'''

def read_crystal_cube_parameters():
	with open('crystal.cube') as crf:
		head = [next(crf) for x in range(6)]
		start = [float(x) for x in head[2].split()[1:]]
		vectors = [0.,0.,0.]
		steps = [0,0,0]
		for i in range(3):
			vector = [float(x) for x in head[i + 3].split()[1:]]
			steps[i] = float(head[i + 3].split()[0])
			for j in range(3):
				if i == j:
					vectors[i] = vector[i]
				else:
					assert vector[j] == 0, 'Only for orthogonal cells, sorry'
		return start, vectors, steps

def create_cube_input(start, vectors, steps, coords, r):
	x0, y0, z0 = [int((coords[i] - r) / vectors[i]) * vectors[i] for i in range(3)]
	x1, y1, z1 = [int((coords[i] + r) / vectors[i]) * vectors[i] for i in range(3)]

	dx, dy, dz = vectors
	nx, ny, nz = [int(([x1, y1, z1][i]-[x0, y0, z0][i])/vectors[i]) for i in range(3)]
	res = cube_template % (x0, y0, z0, nx, ny, nz, dx, dy, dz)
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
				elif line.startswith('$density-cube'):
					skipping = True
				if not skipping:
					outp.write(line)
			
		
r = float(sys.argv[1])
			
start, vectors, steps = read_crystal_cube_parameters()
coords =read_center_coords()

cube_input = create_cube_input(start, vectors, steps, coords, r)

rewrite_control(cube_input)


















