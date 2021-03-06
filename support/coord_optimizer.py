#!/usr/bin/python 
import os, sys, threading, time

def read_embed_positions():
	empos = []
	with open("coord") as c:
		n = 0
		for l in c:
			ls = l.split()
			if '$' not in l and ls:
				if ls[-1] == 'zz':
					empos += [n]
				n += 1
	return n, empos	

def read_grad_from_control(num, empos):
	rewrite = ''
	grad = 0
	with open("control") as c:
		lines = c.read().splitlines()
		n = 0
		while n < len(lines):
			if '$grad' in lines[n]:
				break
			rewrite += lines[n] + '\n'
			n += 1
		assert n < len(lines)
		grad = 0.
		while n < len(lines):
			end = False
			while n < len(lines):
				if lines[n][0] == '#':
					rewrite += lines[n] + '\n'
					n += 1
					continue
				if 'cycle' in lines[n]:
					break
				if '$' in lines[n] and '$grad' not in lines[n]:
					end = True
					break
				rewrite += lines[n] + '\n'
				n += 1
			if n == len(lines) or end:
				break
			grad = 0.
			rawgrads = []
			n0 = n
			for i in range(num):
				rewrite += lines[n] + '\n'
				n = n0 + 1 +  i
			for i in range(num):
				rewrite += lines[n] + '\n'
				n = n0 + num + 1 +  i
				igrads = [float(l.replace("D", "E")) for l in lines[n].split()]
				if i not in empos:
					rawgrads.extend(igrads)
				else:
					lines[n] = '0.0 0.0 0.0'
			grad = (reduce(lambda s, i: s + i**2, rawgrads, 0.))**0.5
		while n < len(lines):
				rewrite += lines[n] + '\n'
				n += 1
	with open("control", "w") as c1:
		c1.write(rewrite)
	return grad

def read_coord(filename):
	with open(filename) as inp:
		coords = []
		lines = inp.read().splitlines()
		assert lines[0] == "$coord"
		assert lines[-1] == "$end"
		for line in lines[1:-1]:
			ls = line.strip().split()
			name = ls[3] if len(ls) == 4 else 'zz'
			x,y,z = [float(i) for i in ls[0:3]]
			coords.append([name, [x,y,z]])
		return coords
		
def write_coords(coords, filename):
	with open(filename, 'w') as outp:
		outp.write("$coord\n")
		for name, [x,y,z] in coords:
			outp.write(" %16.10f %16.10f% 16.10f  %s\n" % (x,y,z,name))
		outp.write("$end")

def write_result(grad, damp):
	with open("coord.log", "a") as log:
		log.write("grad = %13.8f | damp = %5.3f\n" % (grad, damp))

def step(startcoords, lastcoords, num, empos, lastgrad, damp):
		assert os.system('/home/demidov/turbo/bin/amd64/part_relax > log_relax.log') == 0
		time.sleep(2)
		relaxcoords = read_coord('coord')
			
		ncoords = []
		for i, ([sname, sc], [lname, lc], [rname, rc]) in enumerate(zip(startcoords, lastcoords, relaxcoords)):
			nc = None
			if i in empos:
				nc = sc
			else:
				nc = []
				for lk, rk in zip(lc, rc):
					nc.append(lk*(1-damp) + rk*damp)
			ncoords.append([sname, nc])
		write_coords(ncoords, 'coord')
		assert os.system('/home/demidov/turbo/bin/amd64/part_dscf > log_dscf.log') == 0
		time.sleep(2)
		assert os.system('/home/demidov/turbo/bin/amd64/part_grad > log_grad.log') == 0
		time.sleep(2)
		grad = read_grad_from_control(num, empos)
		write_coords(ncoords, 'coord.restart')
		return grad, ncoords

eps = float(sys.argv[1])

coords = None
if os.path.exists('coord.restart'):
	print 'restarting'
	coords = read_coord('coord.restart')
else:
	coords = read_coord('coord.start')
write_coords(coords, 'coord')

num, empos = read_embed_positions()

with open("coord.log", "w") as log:
	log.write("OPTIMIZATION START: eps=%.2e\n" % eps)
assert os.system('/home/demidov/turbo/bin/amd64/part_dscf > log_dscf.log') == 0
time.sleep(2)
assert os.system('/home/demidov/turbo/bin/amd64/part_grad > log_grad.log') == 0
time.sleep(2)
grad = read_grad_from_control(num, empos)
damp = 0.5
write_result(grad, damp)
tcoords = coords

while True:
	grad1, tcoords = step(coords, tcoords, num, empos, grad, damp)
	write_result(grad1, damp)
	if abs(grad1-grad) < eps:
		break
	grad = grad1




























