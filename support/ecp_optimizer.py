#!/usr/bin/python 
import os, sys, math
from scipy.optimize import minimize
from functools import reduce

dscf_command = '/home/maltsev/from_demidov/turbo/bin/amd64/part_dscf>log_dscf.log'
grad_command = '/home/maltsev/from_demidov/turbo/bin/amd64/part_grad>log_grad.log'

last_xe = None
last_xe_step = None

PREC = 2e-8

def check_if_step_not_grad(xe, eps):
	xe = list(xe)
	global last_xe_step
	global last_xe
	if last_xe == None:
		last_xe_step = xe
		last_xe = xe
		return True
	else:
		assert len(xe) == len(last_xe_step)
		n_diff = 0
		s_diff = 0
		for x, lx in zip(xe, last_xe_step):
			if abs(x - lx) >= PREC:
				n_diff += 1
				s_diff += abs(x - lx)
				
		if n_diff == 1 and abs(s_diff - eps) <= PREC:
			last_xe = xe
			return False
		else:
			last_xe_step = xe
			last_xe = xe
			return True

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
	with open("control") as c:
		lines = c.read().splitlines()
		n = 0
		while n < len(lines):
			if '$grad' in lines[n]:
				break
			n += 1
		assert n < len(lines)
		grad = 0.
		while n < len(lines):
			while n < len(lines):
				if lines[n][0] == '#':
					n += 1
					continue
				if 'cycle' in lines[n]:
					break
				n += 1
			if n == len(lines):
				return grad
			grad = 0.
			rawgrads = []
			n0 = n
			for i in range(num):
				n = n0 + num + 1 +  i
				igrads = [float(l.replace("D", "E")) for l in lines[n].split()]
				if i not in empos:
					rawgrads.extend(igrads)
			grad = (reduce(lambda s, i: s + i**2, rawgrads, 0.))**0.5

def read_basis_start(filename):
	with open(filename) as inp:
		basis = inp.read()
		lines = basis.splitlines()
		exps = []
		for line in lines:
			if line and line[-1] == '!':
				exp = math.log(float(line.split()[-2]))
				exps.append(exp)
		return basis, exps

def write_basis(basis, exps, restart):
	with open('basis.restart' if restart else 'basis', 'w') as outp:
		lines = basis.splitlines()
		expindex = 0
		for line in lines:
			if line and line[-1] == '!':
				ls = line.split()
				ls[-2] = math.exp(exps[expindex])
				expindex += 1
				fmt = "   %f   %d   %f  !" if restart else "   %f   %d   %f"
				line = fmt % (float(ls[0]), int(ls[1]), float(ls[2]))
			outp.write(line + "\n")

def write_result(basis, exps, grad, truestep):
	with open("basis.log", "a") as logf:
		logf.write("********** VALUE = %12.10f | TYPE = %s **********\n" % (grad, "MAIN" if truestep else "GRAD"))
		for e in exps:
			logf.write('  log ={:12.7f} value ={:12.7f}\n'.format(e, math.exp(e)))

def calculate(n, empos, basis, exps, eps):
		write_basis(basis, exps, False)
		assert os.system(dscf_command) == 0
		assert os.system(grad_command) == 0
		grad = read_grad_from_control(n, empos)
		t = check_if_step_not_grad(exps, eps)
		write_result(basis, exps, grad, t)
		if t:
			write_basis(basis, exps, True)
		return grad

def optimize_ecp(eps, ftol, maxit):
	filename = 'basis.restart' if os.path.exists('basis.restart') else 'basis.start'
#	filename = 'basis.start'
	n, empos = read_embed_positions()
	basis, exps = read_basis_start(filename)
	target = lambda x: calculate(n, empos, basis, x, eps)
	bounds = [(math.log(1./15),math.log(15))] * len(exps)
	res = minimize(target, exps, bounds=bounds, method='SLSQP', options={'eps':eps, 'ftol':ftol, 'maxiter':maxit})

if len(sys.argv) == 3:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])
	with open("basis.log", "w") as bb:
		bb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))
	optimize_ecp(eps, ftol, 9999)

else:
	print('Usage: ecp_optimizer.py [STEP] [PRECISION]\nExample: charge_optimizer.py 1e-2 1e-5')




























