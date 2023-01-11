#!/usr/bin/env python
import sys, os, re, math, random
import scipy as sc
import subprocess as sp 
from scipy.optimize import minimize
filename = sys.argv[1]
start_filename = filename + '.start'
restart_filename = filename + '.restart'
log_filename = filename + '.ecplog'
d12_filename = filename + '.d12'
out_filename = filename + '.out'
crystall_command = "{crystal}/runcry17 {filename} {filename}".format(crystal=os.environ['CRYSTAL'],filename=filename)

step = float(sys.argv[2])
prec = float(sys.argv[3])

class GroupData:
	def __init__(self, name):
		self.name = name
		self.effective_value = 1
		self.value_is_2 = False

class EcpManager:
	def __init__(self):
		self.groups = []
		self.rawfile = ''

def read_grad(out_filename=out_filename):
	with open(out_filename) as f_in:
		state=0
		for line in f_in:
			if 'CRYSTAL17' in line:
				state = 1
				continue
			elif 'RMS GRADIENT' in line and state==1:
				return float(line.split()[2]) 
		if state == 0:
			raise ValueError('{} is not crystall09 out file'.format(out_filename))
		else:
			raise ValueError('{} is invalid crystall09 out file'.format(out_filename))
		
def read_start():
	filename = restart_filename if os.path.exists(restart_filename) else start_filename
	with open(filename) as inpf:
		e = EcpManager()
		for line in inpf.readlines():
			if '!' in line:
				group_name = line.split('!')[1].strip()
				group = None
				group_templist = [g for g in e.groups if g.name == group_name]
				if group_templist:
					group = group_templist[0]
				else:
					group = GroupData(group_name)
					e.groups.append(group)
				
				a, k, nr = [float(x) for x in line.split('!')[0].split()]
				nr = int(nr)
				if nr == -2:
					group.effective_value = a
					group.value_is_2 = True
				elif not group.value_is_2:
					group.effective_value = a
				

			e.rawfile += line
	return e

def write_current(e, values, restart, noguess):
	filename = restart_filename if restart else d12_filename
	with open(filename, 'w') as outpf:
		for line in e.rawfile.splitlines():
			if noguess and 'GUESSP' in line:
				continue
			if '!' in line:
				group_name = line.split('!')[1].strip()
				group = [g for g in e.groups if g.name == group_name][0]
				group_index = e.groups.index(group)
				value = values[group_index]
				a, k, nr = [float(x) for x in line.split('!')[0].split()]
				nr = int(nr)
				mult = math.exp(value) / group.effective_value
				if nr == -2:
					a *= mult
				elif nr == -1:
					a *= mult
					k *= mult**0.5
				elif nr == 0:
					a *= mult
					k *= mult
				else:
					assert False, nr
				line = "{:24.16f} {:24.16f} {:3d}".format(a, k, nr)
				if restart:
					line += " !" + group_name
			outpf.write(line + '\n')

last_xs = None
last_xs_step = None

def check_if_step_not_grad(xs, eps):
	PREC = 2e-8
	xs = list(xs)
	global last_xs_step
	global last_xs
	if last_xs == None:
		last_xs_step = xs
		last_xs = xs
		return True
	else:
		assert len(xs) == len(last_xs_step)
		n_diff = 0
		s_diff = 0
		for x, lx in zip(xs, last_xs_step):
			if abs(x - lx) >= PREC:
				n_diff += 1
				s_diff += abs(x - lx)
				
		if n_diff == 1 and abs(s_diff - eps) <= PREC:
			last_xs = xs
			return False
		else:
			last_xs_step = xs
			last_xs = xs
			return True
			
def write_log_before(t, values):
	with open(log_filename, 'a') as logf:
		logf.write("TYPE = {}\n".format("MAIN" if t else "GRAD"))
		for v in values:
			logf.write('  log ={:12.7f} value ={:12.7f}\n'.format(v, math.exp(v)))

def write_log_after(grad):
	with open(log_filename, 'a') as logf:
		logf.write("  VALUE = {:20.15f}\n".format(grad))
		
def write_log_error(err):
	with open(log_filename, 'a') as logf:
		logf.write("  %s\n" % err)
		
	
def do_step(xs, e, crystall_command=crystall_command):
	write_current(e, xs, False, False)
	t = check_if_step_not_grad(xs, step)
	write_log_before(t, xs)
	ret = sp.call(crystall_command, shell=True)
	print(ret)
#	with open(out_filename) as outf:
#		if 'SIGTERM' in outf.read():
#			write_log_error("ERROR, restarting without GUESSP...")
#			write_current(e, xs, False, True)
#			sp.call(crystall_command, shell=True)
#			with open(out_filename) as outf1:
#				if 'SIGTERM' in outf1.read():
#					write_log_error("...does not help; exiting :(")
#					assert False
	grad = read_grad()
	if t:
		write_current(e, xs, True, False)
	write_log_after(grad)
	return grad

e = read_start()
values = [math.log(g.effective_value) for g in e.groups]
bounds = [(math.log(1./15),math.log(15))] * len(values)
target = lambda x: do_step(x, e)

with open(log_filename, 'w') as logf:
	logf.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (step, prec))
res = minimize(target, values, bounds=bounds, method='SLSQP', options={'eps':step, 'ftol':prec, 'maxiter':30})
#res = minimize(target, values, bounds=bounds, method='L-BFGS-B', options={'eps':step, 'gtol':prec, 'maxfun':500, 'maxiter':30})
write_current(e, res.x, False, False)
write_current(e, res.x, True, False)
