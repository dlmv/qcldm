
import math, sys, os
from scipy.optimize import minimize

last_xe = None
last_xe_step = None
PREC = 2e-8

shells = 'spdfgh'

def read_ini():
	with open('optbas.ini') as ini:
		optnums = {}
		lines = ini.read().splitlines()
		name = lines[0]
		for line in lines[1:]:
			if line:
				ls = line.split()
				optnums[shells.index(ls[0])] = int(ls[1])
		return name, optnums

class BasisPrimitive:
	def __init__(self, a, c):
		self.a = a
		self.c = c

class BasisFunction:
	def __init__(self, primitives):
		self.primitives = primitives

class Basis:
	def __init__(self, name, comment, data):
		self.data = data
		self.name = name
		self.comment = comment
	
	def tostring(self):
		res = self.name + "\n" + self.comment + "\n*\n"
		for l in sorted(self.data.keys()):
			for f in self.data[l]:
				res += " %3d   %s\n" % (len(f.primitives), shells[l])
				for p in f.primitives:
					res += " %16.10f   %16.10f\n" % (p.a, p.c)
		return res

	def opt_params(self, optnums):
		xs = []
		for l in sorted(optnums.keys()):
			assert len(self.data[l]) >= optnums[l]
			for n in range(optnums[l]):
				assert len(self.data[l][n].primitives) == 1
				xs.append(math.log(self.data[l][n].primitives[0].a))
		return xs

	def update_params(self, optnums, xs):
		i = 0
		for l in sorted(optnums.keys()):
			assert len(self.data[l]) >= optnums[l]
			for n in range(optnums[l]):
				assert len(self.data[l][n].primitives) == 1
				self.data[l][n].primitives[0].a = math.exp(xs[i])
				i += 1
		return xs

	def build_constraints(self, optnums):
		start = 0
		cons = []
		for l in sorted(optnums.keys()):
			if l == 0:
				d = 1.1 
			elif l == 1:
				d = 1.15
			else:
				d = 1.1
			for n in range(optnums[l] - 1):
				cons.append({'type': 'ineq', 'fun' : lambda x,i=n,start=start,d=d: x[start+i] - x[start+i+1] - math.log(d)})
			if len(self.data[l]) > optnums[l] and optnums[l] > 0:
				cons.append({'type': 'ineq', 'fun' : lambda x,i=optnums[l]-1,start=start,d=d,l=l: x[start+i] - math.log(self.data[l][i+1].primitives[0].a) - math.log(d)})
			
			if l == 0:
				d = 1.3**2
				for n in range(optnums[l] - 2):
					cons.append({'type': 'ineq', 'fun' : lambda x,i=n,start=start,d=d: x[start+i] - x[start+i+2] - math.log(d)})
				if len(self.data[l]) > optnums[l] and optnums[l] > 1:
					cons.append({'type': 'ineq', 'fun' : lambda x,i=optnums[l]-2,start=start,d=d,l=l: x[start+i] - math.log(self.data[l][i+2].primitives[0].a) - math.log(d)})
				if len(self.data[l]) > optnums[l] + 1 and optnums[l] > 0:
					cons.append({'type': 'ineq', 'fun' : lambda x,i=optnums[l]-1,start=start,d=d,l=l: x[start+i] - math.log(self.data[l][i+2].primitives[0].a) - math.log(d)})
			start += optnums[l]
		return cons
				
		

def basis_from_string(s):
	lines = s.splitlines()
	name = lines[0]
	comment = lines[1]
	data = {}
	i = 3
	while i < len(lines):
		nc, l = lines[i].split()
		nc, l = int(nc), shells.index(l)
		prims = []
		for i in range(i + 1, i + 1 + nc):
			a,c = lines[i].split()
			prims.append(BasisPrimitive(float(a), float(c)))
		if l not in data.keys():
			data[l] = []
		data[l].append(BasisFunction(prims))
		i += 1
	return Basis(name, comment, data)

class Ecp:	#todo
	def __init__(self, name, raw_data):
		self.name = name
		self.raw_data = raw_data
		
	def tostring(self):
		res = self.name + "\n"
		res += self.raw_data + "\n"
		return res

def ecp_from_string(s):
	lines = s.splitlines()
	name = lines[0]
	raw_data = "\n".join(lines[1:])
	return Ecp(name, raw_data)	

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

def read_energy_from_control():
	with open('control') as control:
		read_energy = False
		last_energy = 0.
		for line in control.xreadlines():
			if '$energy' in line:
				read_energy = True
			elif read_energy and '$' in line:
				read_energy = False
				return last_energy
			elif read_energy:
				last_energy = float(line.split()[1])

def read_start_basis():
	basises = []
	ecps = []
	filename = "basis.start"
	if os.path.exists("basis.restart"):
		filename = "basis.restart"
		print "restarting"
	with open(filename) as basis:
		lines = basis.readlines()
		i = 0
		while "$basis" not in lines[i]:
			i += 1
		i += 1
		while '*' in lines[i] and '$ecp' not in lines[i+1]:
			ast_num = 0
			bas_s = ''
			i += 1
			while True:
				if '*' in lines[i]:
					ast_num += 1
					if ast_num == 2:
						break
				bas_s += lines[i]
				i += 1
			basises.append(basis_from_string(bas_s[:-1]))
		while "$ecp" not in lines[i]:
			i += 1
		i += 1
		while '*' in lines[i] and '$end' not in lines[i+1]:
			ast_num = 0
			ecp_s = ''
			i += 1
			while True:
				if '*' in lines[i]:
					ast_num += 1
					if ast_num == 3:
						break
				ecp_s += lines[i]
				i += 1
			ecps.append(ecp_from_string(ecp_s[:-1]))
		return basises, ecps

def write_basis(basises, ecps, filename='basis'):
	with open(filename, 'w') as basis:
		basis.write("$basis\n")
		for b in basises:
			basis.write("*\n")
			basis.write(b.tostring())
		basis.write("*\n")
		basis.write("$ecp\n")
		for e in ecps:
			basis.write("*\n")
			basis.write(e.tostring())
		basis.write("*\n")
		basis.write("$end")
			
def write_result(e, x, t):
	with open('basis.log', 'a') as logf:
		logf.write("TYPE = {}   VALUE = {:16.10f}\n".format("MAIN" if t else "GRAD", e))
		for v in x:
			logf.write('  log ={:12.7f} value ={:12.7f}\n'.format(v, math.exp(v)))

def calculate(x, optnums, basis, basises, ecps, eps):
	basis.update_params(optnums, x)
	write_basis(basises, ecps)
	
	assert os.system('dscf > log_dscf.log') == 0
	max_dscf = 2
	while True:
		with open('log_dscf.log') as log:
			if 'ENERGY CONVERGED' in log.read():
				break
		assert max_dscf > 0
		assert os.system('dscf > log_dscf.log') == 0
		max_dscf -= 1
	os.system('cp log_dscf.log log_dscf.bak')
	e = read_energy_from_control()

#	e = sum(map(lambda a:a,x))
	
	t = check_if_step_not_grad(x, eps)
	write_result(e, x, t)
	if t:
		write_basis(basises, ecps, 'basis.restart')
	return e

def optimize_basis(eps, ftol, maxit):
	name, optnums = read_ini()
	basises, ecps = read_start_basis()
	basis = next(x for x in basises if x.name == name)
	x0 = basis.opt_params(optnums)
	cons = basis.build_constraints(optnums)

	f = lambda x: calculate(x, optnums, basis, basises, ecps, eps)
	bounds = [(math.log(1./200),math.log(200))] * len(x0)
	res = minimize(f, x0, args=(), method='SLSQP', jac=None, bounds=bounds, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
	




if len(sys.argv) == 3:
	eps = float(sys.argv[1])
	ftol = float(sys.argv[2])
	with open("basis.log", "w") as emb:
		emb.write("OPTIMIZATION START: eps=%.2e, ftol=%.2e\n" % (eps, ftol))
	optimize_basis(eps, ftol, 9999)
	
	



















		
