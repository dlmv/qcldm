#!/usr/bin/python 
import os, sys, threading, time, math
from functools import reduce

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

def write_if_best(grad, coords):
	last_grad = 1e10
	filename = 'coord.best'
	if os.path.exists(filename):
		with open(filename) as cb:
			last_grad = float(cb.readlines()[0])
	if grad < last_grad:
		with open(filename, 'w') as outp:
			outp.write("%13.8f\n" % grad)
			outp.write("$coord\n")
			for name, [x,y,z] in coords:
				outp.write(" %16.10f %16.10f% 16.10f  %s\n" % (x,y,z,name))
			outp.write("$end")

def write_result(grad, scale):
	with open("coord.log", "a") as log:
		log.write("grad = %13.8f | scale = %5.3f\n" % (grad, scale))

def step(startcoords, lastcoords, num, empos, lastgrad, scale):
		assert os.system('%s/bin/amd64/part_relax > log_relax.log' % os.environ['TURBO']) == 0
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
					nc.append(lk*(1-scale) + rk*scale)
			ncoords.append([sname, nc])
		write_coords(ncoords, 'coord')
		assert os.system('%s/bin/amd64/part_dscf > log_dscf.log'  % os.environ['TURBO']) == 0
		time.sleep(2)
		assert os.system('%s/bin/amd64/part_grad > log_grad.log' % os.environ['TURBO']) == 0
		time.sleep(2)
		grad = read_grad_from_control(num, empos)
		write_coords(ncoords, 'coord.restart')
		write_if_best(grad, ncoords)
		return grad, ncoords

eps = float(sys.argv[1])
scale = float(sys.argv[2]) if len(sys.argv) > 2 else 1


coords = None
if os.path.exists('coord.restart'):
	print('restarting')
	coords = read_coord('coord.restart')
else:
	if not os.path.exists('coord.start'):
		os.system('cp coord coord.start')
	coords = read_coord('coord.start')
write_coords(coords, 'coord')

num, empos = read_embed_positions()

with open("coord.log", "w") as log:
	log.write("OPTIMIZATION START: eps=%.2e\n" % eps)
assert os.system('%s/bin/amd64/part_dscf > log_dscf.log' % os.environ['TURBO']) == 0
time.sleep(2)
assert os.system('%s/bin/amd64/part_grad > log_grad.log' % os.environ['TURBO']) == 0
time.sleep(2)
grad = read_grad_from_control(num, empos)
write_result(grad, scale)
tcoords = coords

order = 0
order_count = 0
same_order_limit = 15

while True:
	grad, tcoords = step(coords, tcoords, num, empos, grad, scale)
	write_result(grad, scale)
	if grad < eps:
		break
	grad_order = math.floor(math.log(grad, 10))
	if grad_order == order:
		order_count += 1
		if order_count > same_order_limit:
			break
	else:
		order = grad_order
		order_count = 0




























