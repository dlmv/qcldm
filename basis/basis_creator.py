import os, sys, re, logging, io, shutil, time
import json, zipfile, csv
import math, scipy, numpy
import pyscf

from urllib.request import urlopen
from urllib.parse import quote
from operator import itemgetter
from functools import reduce



        #=======================================#
#=======#                 UTILS                 #=======#
        #=======================================#

def add_coloring_to_emit_ansi(fn):
    # add methods we need to the class
    def new(*args):
        levelno = args[1].levelno
        if(levelno >= 50):
            color = '\x1b[31m' # red
        elif(levelno >= 40):
            color = '\x1b[31m' # red
        elif(levelno >= 30):
            color = '\x1b[33m' # yellow
        elif(levelno >= 20):
            color = '\x1b[32m' # green 
        elif(levelno >= 10):
            color = '\x1b[35m' # pink
        else:
            color = '\x1b[0m' # normal
        args[1].msg = color + args[1].msg + '\x1b[0m'  # normal
        #print "after"
        return fn(*args)
    return new
def colorize():
	logging.StreamHandler.emit = add_coloring_to_emit_ansi(logging.StreamHandler.emit)

logging.basicConfig(format = '%(asctime)s %(levelname).1s %(message)s', level = logging.DEBUG,
                              datefmt='%Y-%m-%d %H:%M:%S')
colorize()

def halt(condition, message):
	if condition:
		logging.error(message)
		assert False

def yesno(message):
	result = False
	while True:
		logging.info(' %s (y/n)' % message)
		s = input().lower().strip()
		if s == 'y':
			result = True
			break
		elif s == 'n':
			break
		else:
			logging.warning(' Input not recognized. Try again.')
	return result

def multiselect(message, variants):
	while True:
		logging.info(message)
		for i, c in enumerate(variants):
			logging.info(' %2d. %s' % (i+1, c))
		s = input().lower().strip()
		if re.match(r"^\d+$", s) and int(s) > 0 and int(s) <= len(variants):
			index = int(s) - 1
			return index, variants[index]
			break
		else:
			logging.warning(' Input not recognized, try again:')


class InputVariant:
	def __init__(self, hint, regex, callback):
		self.hint = hint
		self.regex = regex
		self.callback = callback #return True if accepted, return False to continue input cycle
		self.active = True
		self.custom_field = None


def complexinput(message, variants):
	while True:
		logging.info(message)
		for v in variants:
			if v.active:
				logging.info("*%s" % v.hint)
		s = input().lower().strip()
		variant = None
		for v in variants:
			if v.active:
				if re.match(v.regex, s):
					variant = v
					break
		if variant:
			if variant.callback(s, variant):
				return
		else:
			logging.warning(' Input not recognized, try again:')

def print_title(s):
	logging.info('')
	logging.info('=' * (len(s) + 4))
	logging.info('| %s |' % s)
	logging.info('=' * (len(s) + 4))

def print_subtitle(s):
	logging.info('')
	logging.info('%s' % s)
	logging.info('-' * (len(s)))

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


        #=======================================#
#=======#               ELEMENTS                #=======#
        #=======================================#

SHELLS = 'spdfghi'

ELEMENTS = 'n', \
           'H', 'He', \
           'Li', 'Be',  'B',  'C',  'N',  'O',  'F', 'Ne', \
           'Na', 'Mg', 'Al', 'Si',  'P',  'S', 'Cl', 'Ar', \
           'K',  'Ca', 'Sc', 'Ti',  'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', \
           'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', \
           'Rb', 'Sr',  'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', \
           'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te',  'I', 'Xe', \
           'Cs', 'Ba', 'La', \
           'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', \
                             'Hf', 'Ta',  'W', 'Re', 'Os', 'Ir', 'Pt', \
           'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', \
           'Fr', 'Ra', 'Ac', \
           'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', \
                             'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', \
           'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og'

def get_element(name):
	name = name[0].upper() + name[1:].lower()
	if name in ELEMENTS:
		n = ELEMENTS.index(name)
		return name, n
	else:
		halt(True, "...ELEMENT NOT FOUND: '%s'" % name)

def parse_element_string(s):
	n = -1
	logging.debug("Input string: '%s'" % s)
	if re.match(r'[a-z]{1,2}$', s, re.I):
		logging.debug("...recognized as element name: '%s'" % s)
		s, n = get_element(s)
	elif re.match(r'\d{1,3}$', s, re.I):
		n = int(s)
		logging.debug("...recognized as element number: '%d'" % n)
		halt((n >= len(ELEMENTS)) or (n <= 0), "...ELEMENT IS NOT DEFINED FOR NUMBER %d" % n)
	else:
		halt(True, "...INPUT NOT RECOGNIZED: '%s'" % s)
	logging.debug("Recognized element: %s/%d" % (ELEMENTS[n], n))
	return n

def aufbau_sequence():
	nl = 1
	while True:
		for l in reversed(list(range(nl))):
			n = nl - l
			if n > l:
				yield n, l
		nl += 1

def shell_sequence():
	n = 1
	while True:
		for l in range(n):
			yield n, l
		n += 1

def estimate_configuration(atom_number):
	configuration = {}
	i = 0
	for n, l in aufbau_sequence():
		shell_occ = 4*l + 2
		if i + shell_occ < atom_number:
			configuration[(n,l)] = shell_occ
			i += shell_occ
		else:
			configuration[(n,l)] = atom_number - i
			break
	return configuration

def estimate_valence_configuration(atom_number):
	configuration = estimate_configuration(atom_number)
	nmax = 0
	valence_configuration = []
	for n, l in aufbau_sequence():
		if not (n, l) in configuration.keys():
			break
		if n > nmax:
			valence_configuration = []
		nmax = max(n, nmax)
		if l == 0 or configuration[(n,l)] != 4*l + 2:
			valence_configuration.append([n, l, configuration[(n,l)]])
		if l == 1 and configuration[(n,l)] == 6:
			valence_configuration = []
	return valence_configuration

def outer_core_ns(atom_number):
	configuration_map = estimate_configuration(atom_number)
	valence_configuration = estimate_valence_configuration(atom_number)
	valence_configuration_map = {}
	for n, l, occ in valence_configuration:
		valence_configuration_map[(n, l)] = occ
	max_core_ns = {}
	for (n, l) in configuration_map.keys():
		if (n, l) in valence_configuration_map:
			continue
		if l not in max_core_ns.keys():
			max_core_ns[l] = n
		else:
			max_core_ns[l] = max(max_core_ns[l], n)
	return max_core_ns

def estimate_ion_valence_shells(atom_number, charge=0):
	valence_configuration = estimate_valence_configuration(atom_number)
	if charge:
		valence_configuration_map = {}
		for n, l, occ in valence_configuration:
			valence_configuration_map[(n, l)] = occ
		keys = valence_configuration_map.keys()
		for e in range(charge):
			for n, l in sorted(keys, key=itemgetter(0,1), reverse=True):
				if valence_configuration_map[(n, l)] > 0:
					valence_configuration_map[(n, l)] -= 1
					break
		for conf_part in valence_configuration:
			n, l, occ = conf_part
			conf_part[2] = valence_configuration_map[(n, l)]
	spin = 0
	valence = 0
	valence_string = ''
	for n, l, occ in valence_configuration:
		if occ:
			s = min(occ, 4*l + 2 - occ)
			spin += s
			valence += occ
			valence_string += '.%d%s%d' % (n, SHELLS[l], occ)
	core_num = atom_number - valence - charge
	valence_string = '[%s]' % ELEMENTS[core_num] + valence_string
	return valence, valence_string, spin

def estimate_ecp_configuration(atom_number, ncore):
	configuration = estimate_configuration(atom_number)
	electrons_processed = 0
	ncore_found = ncore == 0
	maxn = 0
	ecp_configuration = []
	failed = False
	for (n, l) in configuration.keys():
		maxn = max(maxn, n)
	for n, l in shell_sequence():
		if n > maxn:
			break
		if (n, l) in configuration.keys():
			if ncore_found:
				ecp_configuration.append([n, l, configuration[(n,l)]])
			else:
				electrons_processed += configuration[(n, l)]
				if electrons_processed == ncore:
					ncore_found = True
				elif electrons_processed > ncore:
					failed = True
					logging.debug("Shell sequence failed, switching to Aufbau")
					break
	if failed:
		ecp_configuration = []
		electrons_processed = 0
		for n, l in aufbau_sequence():
			if n > maxn:
				break
			if (n, l) in configuration.keys():
				if ncore_found:
					ecp_configuration.append([n, l, configuration[(n,l)]])
				else:
					electrons_processed += configuration[(n, l)]
					if electrons_processed == ncore:
						ncore_found = True
					elif electrons_processed > ncore:
						logging.warning("Unexpected NCore: %d" % ncore)
						return None
	print(ecp_configuration)
	return ecp_configuration

def estimate_maxl(atom_number, ncore):
	ecp_configuration = estimate_ecp_configuration(atom_number, ncore)
	if ecp_configuration is None:
		logging.warning("Unexpected NCore: %d; geting maxL from full electron configuration" % ncore)
		return estimate_maxl(atom_number, 0)
	maxl = 0
	for (n, l, occ) in ecp_configuration:
		maxl = max(maxl, l)
	return maxl

def lower_ecp_ns(atom_number, ncore):
	ns = {}
	ecp_configuration = estimate_ecp_configuration(atom_number, ncore)
	if ecp_configuration is None:
		return None
	for n, l, occ in ecp_configuration:
		if l in ns.keys():
			ns[l] = min(ns[l], n)
		else:
			ns[l] = n
	return ns





        #=======================================#
#=======#               GAUSSIANS               #=======#
        #=======================================#

def dfac(n):
	return 1 if n < 2 else reduce(lambda x,y: y*x, list(range(n,1,-2)))

def gauss_norm(a, l):
	return (2**(2*l+3.5)  / dfac(2*l+1) / math.pi**0.5)**0.5 * a**((2.*l+3)/4)


class GaussFunction:
	def __init__(self, a, l):
		self.a = a
		self.l = l

	def __call__(self, r):
		return (r**self.l) * math.exp(-self.a * r**2)


class GaussFunctionNormed(GaussFunction):
	def __init__(self, a, l):
		GaussFunction.__init__(self, a,l)
		self.norm = gauss_norm(a, l)

	def __call__(self, r):
		return GaussFunction.__call__(self, r) * self.norm


def overlap(f1, f2):
	if f1.l != f2.l:
		return 0
	l = f1.l
	ff = GaussFunctionNormed((f1.a + f2.a) / 2, l)
	return (f1.norm * f2.norm) / ff.norm**2

class GaussFunctionContracted:
	def __init__(self, primitives):
		self.primitives = primitives

	def __call__(self, r):
		res = 0
		for c, f in self.primitives:
			res += f(r) * c
		return res

	def get_maximum_a(self):
		r = self.get_maximum_r()
		a = (self.primitives[0][1].l + 1) / r**2 / 2
		return a

	def get_maximum_r(self):
		loggrid = list(numpy.logspace(math.log(0.0001), math.log(10), 500, True, math.exp(1)))
		max_v = 0
		max_r = 0
		for r in loggrid:
			v = abs(self.__call__(r) * r)
			if v > max_v:
				max_v = v
				max_r = r
		return max_r

	def self_overlap(self):
		olp_sum = 0
		for c1, f1 in self.primitives:
			for c2, f2 in self.primitives:
				olp = c1 * c2 * overlap(f1, f2)
				olp_sum += olp
		return olp_sum

	def normalize(self):
		norm = self.self_overlap()
		for cf in self.primitives:
			cf[0] /= norm ** 0.5


	def get_coeff_variables(self):
		maxc = 0
		maxi = -1
		for i, (c,_) in enumerate(self.primitives):
			if abs(c) > abs(maxc):
				maxc = c
				maxi = i
		coeffs = []
		for i, (c,_) in enumerate(self.primitives):
			if i != maxi:
				coeffs.append(c)
		return coeffs, maxi

	def set_coeff_variables(self, coeffs, maxi):
		i = 0
		for c in coeffs:
			if i == maxi:
				i += 1
			self.primitives[i][0] = c
			i += 1

	def get_exp_variables(self):
		exps = []
		for _,f in self.primitives:
			exps.append(math.log(f.a))
		return exps

	def set_exp_variables(self, exps):
		for (_,f), exp in zip(self.primitives, exps):
			f.a = math.exp(exp)

class Projector:
	def __init__(self, core, u):
		self.core = core
		self.u = u

        #=======================================#
#=======#              ECP & BASIS              #=======#
        #=======================================#

class GRECP:
	def __init__(self, n, nval, local, semilocals, spinorbit, projectors):
		self.n = n
		self.nval = nval
		self.ncore = n - nval
		self.local = local #cg
		self.semilocals = semilocals #{l: cg}
		self.spinorbit = spinorbit #{l: cg}
		self.projectors = projectors #{l: [p]}

class Basis:
	def __init__(self, n, functions):
		self.n = n
		self.functions = functions #{l: [cg]}

	def is_uncontracted(self):
		for l in self.functions.keys():
			for cgf in self.functions[l]:
				if len(cgf.primitives) > 1:
					return False
		return True

	def uncontract(self):
		new_functions = {}
		for l in self.functions.keys():
			exponents = set()
			for cgf in self.functions[l]:
				for _, gf in cgf.primitives:
					exponents.add(gf.a)
			new_functions[l] = []
			for e in sorted(list(exponents), reverse=True):
				cgf = GaussFunctionContracted([[1, GaussFunctionNormed(e, l)]])
				new_functions[l].append(cgf)
		self.functions = new_functions

	def maxl(self):
		maxl = 0
		for l in self.functions.keys():
			maxl = max(maxl, l)
		return maxl

	def cut_polarization(self, maxl):
		for l in list(self.functions.keys()):
			if l > maxl:
				del self.functions[l]

def estimate_core_bounds(basis, ncore): #all-electron basis should be used!
	ns = lower_ecp_ns(basis.n, ncore)
	bounds = {}
	for l in sorted(ns.keys()):
		ovcf = basis.functions[l][ns[l] - l - 1]
		bounds[l] = ovcf.get_maximum_r()
	return bounds



        #=======================================#
#=======#             GRECP READING             #=======#
        #=======================================#

def read_table_part(lines, is_ecp, n):
	header = lines[n]
	prim_count = int(header.split()[0])
	proj_count = int(header.split()[1])
	headers = []
	for m in re.finditer(r'\d([SPDFGHI])[^\s]+', header):
		title = m.group()
		l = SHELLS.index(m.group(1).lower())
		span = m.span()
		headers.append((title, l))
	_,L = headers[0]
	for _,l in headers:
		halt(l != L, "DIFFERENT L IN SAME TABLE PART:\n'%s'" % header)
	funcs = {}
	if is_ecp:
		halt(proj_count + 1 * (1 if L == 0 else 2) != len(headers), "WRONG NUMBER OF ECP TITLES:\n'%s'" % header)
	else:
		halt(proj_count != len(headers), "WRONG NUMBER OF FUNC TITLES:\n'%s'" % header)
	for title,_ in headers:
		funcs[title] = GaussFunctionContracted([])
	for n in range(n + 1, n + 1 + prim_count):
		line = lines[n]
		ls = [float(x) for x in line.split()]
		l = L
		if is_ecp:
			l = int(ls[0])
			ls = ls[1:]
		a = ls[0]
		ls = ls[1:]
		if is_ecp and L == 0:
			ls = [ls[0]] + ls[2:]
		expected_length = len(headers)
		halt(len(ls) != expected_length, "WRONG LINE FORMAT:\n'%s'" % line)
		gf = None
		for c, (title,_) in zip(ls, headers):
			if c == 0:
				continue
			gf = GaussFunction(a, l) if is_ecp else GaussFunctionNormed(a, l)
			funcs[title].primitives.append([c, gf])
	return headers, funcs, n + 1

def parse_ecp_mos(mos_format):
	logging.debug("Reading GRECP...")
	lines = mos_format.splitlines()
	m = re.match(r'^([A-Za-z]+)(\d+)[^\W\d]*\s+\*+$', lines[0])
	halt(m is None, "COULD NOT PARSE ECP TITLE: '%s'" % lines[0])
	name = m.group(1)
	nval = int(m.group(2))
	name, nat = get_element(name)
	m = re.match(r'^\s+(\d+(?:\.\d+)?)\s([\d\s]+)$', lines[3])
	halt(m is None, "COULD NOT PARSE ECP HEADER: '%s'" % lines[3])
	nval1f = float(m.group(1))
	nval1 = int(nval1f)
	halt(nval1 != nval1f, "NON-INTEGER NVAL NUMBER: %f" % (nval1f))
	halt(nval != nval1, "DIFFERENT NVAL NUMBERS: %d and %d" % (nval, nval1))
	ls = [int(x) for x in m.group(2).split()]
	nproj = ls[1]
	narep = ls[-3]
	nesop = ls[-2]
	halt(narep != nesop, "DIFFERENT AREP AND ESOP NUMBERS: %d and %d" % (narep, nesop))
	halt(ls[0] != 1, "ATOMS COUNT IS %d INSTEAD OF 1" % ls[0])
	halt(ls[2:-4] != [1] * nproj, "UNEXPECTED CONTRACTION PARAMS: %s INSTEAD OF %s" % (ls[2:-4], [1] * nproj))
	halt(ls[-4] != 1, "CONTRACTION TYPE IS %d INSTEAD OF 1" % ls[0])
	logging.debug("...ECP is for %s (%d/%d)" % (name, nval, nat))
	n = 4
	projector_map = {}
	for i in range(nproj):
		headers, funcs, n = read_table_part(lines, False, n)
		projector_map = {**projector_map, **funcs}
	logging.debug("Projectors: %s" % sorted(projector_map.keys()))
	ecp_map = {}
	for i in range(narep):
		headers, funcs, n = read_table_part(lines, True, n)
		ecp_map = {**ecp_map, **funcs}
	logging.debug("Potentials: %s" % sorted(ecp_map.keys()))
	semilocal_parts = {}
	spinorbit_parts = {}
	for title in ecp_map.keys():
		if m := re.match(r'^\d([SPDFGHI])-(AREP|ESOP)$', title):
			l = SHELLS.index(m.group(1).lower())
			if m.group(2) == 'AREP':
				semilocal_parts[l] = ecp_map[title]
			else:
				spinorbit_parts[l] = ecp_map[title]
		else:
			pass #TODO: projectors
	logging.debug("SemiLocal potentials for Ls: %s" % sorted(semilocal_parts.keys()))
	logging.debug("SpinOrbit potentials for Ls: %s" % sorted(spinorbit_parts.keys()))
	for l in spinorbit_parts.keys():
		for ls in spinorbit_parts[l].primitives:
			ls[0] = ls[0] / (2 * l + 1) * 2
	lrange = list(range(len(semilocal_parts)))
	lrange_so = list(range(1, len(spinorbit_parts) + 1))
	halt(sorted(semilocal_parts.keys()) != lrange, "NON-SEQUENTIAL AREP L RANGE: %s" % sorted(semilocal_parts.keys()))
	halt(sorted(spinorbit_parts.keys()) != lrange_so, "NON-SEQUENTIAL ESOP L RANGE: %s" % sorted(spinorbit_parts.keys()))
	halt(len(lrange_so) + 1 > len(lrange), "ESOP L RANGE EXCEED AREP L RANGE: %s > %s" % (lrange_so, lrange))
	lmax = lrange[-1]
	local_part = semilocal_parts[lmax]
	del semilocal_parts[lmax]
	grecp = GRECP(nat, nval, local_part, semilocal_parts, spinorbit_parts, [])
	return grecp

def read_ecp_mos(filename):
	with open(filename) as f:
		element = None
		return parse_ecp_mos(f.read())

        #=======================================#
#=======#          NWCHEM/PYSCF WRITING         #=======#
        #=======================================#


def ecp_to_nw(ecp):
	res = ''
	res += 'ECP\n'
	res += '%s nelec %d\n' % (ELEMENTS[ecp.n], ecp.n - ecp.nval)
	res += '%s ul\n' % (ELEMENTS[ecp.n])
	for c, gf in ecp.local.primitives:
		res += '%5d %20.8f %20.8f\n' % (gf.l, gf.a, c)
	for l in sorted(ecp.semilocals.keys()):
		res += '%s %s\n' % (ELEMENTS[ecp.n], SHELLS[l].upper())
		cgf = ecp.semilocals[l]
		for c, gf in cgf.primitives:
			res += '%5d %20.8f %20.8f\n' % (gf.l, gf.a, c)
	res += 'END\n'
	if ecp.spinorbit:
		res += '\n'
		res += 'SO\n'
		for l in sorted(ecp.spinorbit.keys()):
			res += '%s %s\n' % (ELEMENTS[ecp.n], SHELLS[l].upper())
			cgf = ecp.spinorbit[l]
			for c, gf in cgf.primitives:
				res += '%5d %20.8f %20.8f\n' % (gf.l, gf.a, c)
		res += 'END\n'
	return res

def write_ecp_nw(ecp, filename):
	with open(filename, 'w') as f:
		f.write(ecp_to_nw(ecp))

def ecp_to_pyscf(ecp):
	res = ''
	res += 'ECP\n'
	res += '%s nelec %d\n' % (ELEMENTS[ecp.n], ecp.n - ecp.nval)
	res += '%s ul\n' % (ELEMENTS[ecp.n])
	l = len(ecp.semilocals)
	cgf = ecp.local
	so_cgf = ecp.spinorbit[l] if l in ecp.spinorbit.keys() else None
	if not so_cgf:
		for c, gf in cgf.primitives:
			res += '%5d %20.8f %20.8f\n' % (gf.l, gf.a, c)
	else:
		exponents = []
		for c, gf in cgf.primitives:
			exponents.append((gf.a, gf.l))
		for c, gf in so_cgf.primitives:
			if (gf.a, gf.l) not in exponents:
				exponents.append((gf.a, gf.l))
		for a, l in exponents:
			c_l = 0
			c_s = 0
			for c, gf in cgf.primitives:
				if gf.a == a and gf.l == l:
						c_l = c
			for c, gf in so_cgf.primitives:
				if gf.a == a and gf.l == l:
						c_s = c
			res += '%5d %20.8f %20.8f %20.8f\n' % (l, a, c_l, c_s)
	for l in sorted(ecp.semilocals.keys()):
		res += '%s %s\n' % (ELEMENTS[ecp.n], SHELLS[l].upper())
		cgf = ecp.semilocals[l]
		so_cgf = ecp.spinorbit[l] if l in ecp.spinorbit.keys() else None
		if not so_cgf:
			for c, gf in cgf.primitives:
				res += '%5d %20.8f %20.8f\n' % (gf.l, gf.a, c)
		else:
			exponents = []
			for c, gf in cgf.primitives:
				exponents.append((gf.a, gf.l))
			for c, gf in so_cgf.primitives:
				if (gf.a, gf.l) not in exponents:
					exponents.append((gf.a, gf.l))
			for a, l in exponents:
				c_l = 0
				c_s = 0
				for c, gf in cgf.primitives:
					if gf.a == a and gf.l == l:
							c_l = c
				for c, gf in so_cgf.primitives:
					if gf.a == a and gf.l == l:
							c_s = c
				res += '%5d %20.8f %20.8f %20.8f\n' % (l, a, c_l, c_s)
	res += 'END\n'
	return res

def write_ecp_pyscf(ecp, filename):
	with open(filename, 'w') as f:
		f.write(ecp_to_pyscf(ecp))

def basis_to_nw(basis):
	res = ''
	res += 'BASIS SPHERICAL\n'
	for l in sorted(basis.functions.keys()):
		for cgf in basis.functions[l]:
			res += '%s   %s\n' % (ELEMENTS[basis.n], SHELLS[l].upper())
			for c, gf in cgf.primitives:
				res += '%20.8f %20.8f\n' % (gf.a, c)
	res += 'END\n'
	return res

def write_basis_nw(basis, filename):
	with open(filename, 'w') as f:
		f.write(basis_to_nw(basis))

        #=======================================#
#=======#             CRYSTAL WRITING           #=======#
        #=======================================#


def ecp_and_basis_to_crystal(ecp, basis, occs):
	basis_size = 0
	for l in basis.functions.keys():
		basis_size += len(basis.functions[l])
	res = "2%02d %d\n" % (ecp.n, basis_size)
	res += "INPUT\n"
	ecp_sizes = [0] * 6
	ecp_sizes[0] = len(ecp.local.primitives)
	for l in sorted(ecp.semilocals.keys()):
		if l > 4:
			break
		ecp_sizes[l + 1] = len(ecp.semilocals[l].primitives)
	res += "%d. %s\n" % (ecp.nval, " ".join([str(n) for n in ecp_sizes]))
	for c, f in ecp.local.primitives:
		res += '%20.8f %20.8f %5d\n' % (f.a, c, f.l - 2)
	for l in sorted(ecp.semilocals.keys()):
		for c, f in ecp.semilocals[l].primitives:
			res += '%20.8f %20.8f %5d\n' % (f.a, c, f.l - 2)
	for l in sorted(basis.functions.keys()):
		ltyp = 0 if l == 0 else l + 1
		for i, cgf in enumerate(basis.functions[l]):
			occ = 0
			if l in occs.keys() and i < len(occs[l]):
				occ = occs[l][i]
			res += '%d %d %d %.1f %.1f\n' % (0, ltyp, len(cgf.primitives), occ, 1)
			for c, f in cgf.primitives:
				res += '%20.8f %20.8f\n' % (f.a, c)
	return res

def write_ecp_and_basis_crystal(ecp, basis, occs, filename):
	with open(filename, 'w') as f:
		f.write(ecp_and_basis_to_crystal(ecp, basis, occs))

        #=======================================#
#=======#          NWCHEM/PYSCF READING         #=======#
        #=======================================#

def parse_ecp_nw(s):
	lines = s.splitlines()
	halt(lines[0].strip() != 'ECP', "UNEXPECTED NW ECP HEADER: '%s'" % (lines[0]))
	grecp = None
	atom_name = None
	if m := re.match(r'^([A-Za-z]+)\s+nelec\s+(\d+)$', lines[1]):
		atom_name, atom_number = get_element(m.group(1))
		ncore = int(m.group(2))
		nval = atom_number - ncore
		semilocal_parts = {}
		local_part = GaussFunctionContracted([])
		n = 2
		l = -1
		while True:
			if lines[n].strip() == 'END':
				grecp = GRECP(atom_number, nval, local_part, semilocal_parts, {}, [])
				n += 1
				break
			else:
				if m := re.match(r'^([A-Za-z]+)\s+([SPDFGHI]+|ul)$', lines[n]):
					atom_name1, _ = get_element(m.group(1))
					halt(atom_name1 != atom_name, "DIFFERENT ELEMENTS IN HEADER: %s and %s" % (atom_name, atom_name1))
					if m.group(2) == 'ul':
						l = -1
					else:
						l = SHELLS.index(m.group(2).lower())
						cg = GaussFunctionContracted([])
						semilocal_parts[l] = cg
				else:
					try:
						ls = lines[n].split()
						ll = int(ls[0])
						a = float(ls[1])
						c = float(ls[2])
						gf = GaussFunction(a, ll)
						if l == -1:
							local_part.primitives.append([c, gf])
						else:
							semilocal_parts[l].primitives.append([c, gf])
					except Exception:
						halt(True, "CANNOT PARSE ECP LINE: '%s'" % lines[n])
			n += 1
		halt(n >= len(lines), "UNEXPECTED ECP END")
	else:
		halt(True, "CANNOT PARSE ECP HEADER: '%s'" % lines[1])
	while n < len(lines) and not (lines[n]).strip():
		n += 1
	if n < len(lines) and (lines[n]).strip() == 'SO':
		n += 1
		spinorbit_parts = {}
		while True:
			if lines[n].strip() == 'END':
				grecp.spinorbit = spinorbit_parts
				break
			else:
				if m := re.match(r'^([A-Za-z]+)\s+([SPDFGHI]+)$', lines[n]):
					atom_name1, _ = get_element(m.group(1))
					halt(atom_name1 != atom_name, "DIFFERENT ELEMENTS IN HEADER: %s and %s" % (atom_name, atom_name1))
					l = SHELLS.index(m.group(2).lower())
					cg = GaussFunctionContracted([])
					spinorbit_parts[l] = cg
				else:
					try:
						ls = lines[n].split()
						ll = int(ls[0])
						a = float(ls[1])
						c = float(ls[2])
						gf = GaussFunction(a, ll)
						spinorbit_parts[l].primitives.append([c, gf])
					except Exception:
						halt(True, "CANNOT PARSE ECP LINE: '%s'" % lines[n])
			n += 1
		halt(n >= len(lines), "UNEXPECTED ECP END")
	return grecp

def read_ecp_nw(filename):
	with open(filename) as f:
		element = None
		return parse_ecp_nw(f.read())

def parse_ecp_pyscf(s):
	lines = s.splitlines()
	halt(lines[0].strip() != 'ECP', "UNEXPECTED NW ECP HEADER: '%s'" % (lines[0]))
	grecp = None
	atom_name = None
	if m := re.match(r'^([A-Za-z]+)\s+nelec\s+(\d+)$', lines[1]):
		atom_name, atom_number = get_element(m.group(1))
		ncore = int(m.group(2))
		nval = atom_number - ncore
		semilocal_parts = {}
		spinorbit_parts = {}
		local_part = GaussFunctionContracted([])
		n = 2
		l = -1
		while True:
			if lines[n].strip() == 'END':
				spinorbit_parts[len(semilocal_parts)] = spinorbit_parts[-1]
				del spinorbit_parts[-1]
				grecp = GRECP(atom_number, nval, local_part, semilocal_parts, spinorbit_parts, [])
				break
			else:
				if m := re.match(r'^([A-Za-z]+)\s+([SPDFGHI]+|ul)$', lines[n]):
					atom_name1, _ = get_element(m.group(1))
					halt(atom_name1 != atom_name, "DIFFERENT ELEMENTS IN HEADER: %s and %s" % (atom_name, atom_name1))
					if m.group(2) == 'ul':
						l = -1
					else:
						l = SHELLS.index(m.group(2).lower())
						cg = GaussFunctionContracted([])
						semilocal_parts[l] = cg
				else:
					try:
						ls = lines[n].split()
						ll = int(ls[0])
						a = float(ls[1])
						c = float(ls[2])
						gf = GaussFunction(a, ll)
						if c != 0:
							if l == -1:
								local_part.primitives.append([c, gf])
							else:
								semilocal_parts[l].primitives.append([c, gf])
						if len(ls) > 3:
							c = float(ls[3])
							if c != 0:
								gf = GaussFunction(a, ll)
								if l not in spinorbit_parts.keys():
									spinorbit_parts[l] = GaussFunctionContracted([])
								spinorbit_parts[l].primitives.append([c, gf])
					except Exception:
						halt(True, "CANNOT PARSE ECP LINE: '%s'" % lines[n])
			n += 1
		halt(n >= len(lines), "UNEXPECTED ECP END")
	else:
		halt(True, "CANNOT PARSE ECP HEADER: '%s'" % lines[1])
	return grecp

def read_ecp_pyscf(filename):
	with open(filename) as f:
		element = None
		return parse_ecp_pyscf(f.read())

def parse_basis_nw(s):
	lines = s.splitlines()
	while not lines[0].strip() or lines[0].strip().startswith('#'):
		del lines[0]
	halt(not re.match(r'^BASIS (.+ )?SPHERICAL.*?', lines[0]), "UNEXPECTED NW BASIS HEADER: '%s'" % lines[0])
	functions = {}
	cgs = []
	l = -1
	n = 1
	while lines[n].strip().startswith('#'):
		n += 1
	atom_name = None
	atom_number = None
	while True:
		if lines[n].strip() == 'END':
			halt(l == -1, "EMPTY BASIS")
			halt(len(cgs) == 0, "EMPTY BASIS FUNCTION")
			if l not in functions.keys():
				functions[l] = []
			for cg in cgs:
				functions[l].append(cg)
			return Basis(atom_number, functions)
		else:
			if m := re.match(r'^([A-Za-z]+)\s+([SPDFGHI]+)$', lines[n]):
				atom_name1, atom_number = get_element(m.group(1))
				if atom_name:
					halt(atom_name1 != atom_name, "DIFFERENT ELEMENTS IN HEADER: %s and %s" % (atom_name, atom_name1))
				atom_name = atom_name1
				if l != -1:
					if l not in functions.keys():
						functions[l] = []
					for cg in cgs:
						functions[l].append(cg)
				l = SHELLS.index(m.group(2).lower())
				cgs = []
			else:
				try:
					ls = lines[n].split()
					a = float(ls[0])
					cs = [float(x) for x in ls[1:]]
					if len(cgs) == 0:
						for c in cs:
							cgs.append(GaussFunctionContracted([]))
					halt(len(cgs) != len(cs), "WRONG NUMBER OF COLUMNS: %s" % lines[n])
					for c, cg in zip(cs, cgs):
						gf = GaussFunctionNormed(a, l)
						cg.primitives.append([c, gf])
				except Exception:
					halt(True, "CANNOT PARSE BASIS LINE: '%s'" % lines[n])
		n += 1
		halt(n >= len(lines), "UNEXPECTED BASIS END")

def read_basis_nw(filename):
	with open(filename) as f:
		element = None
		return parse_basis_nw(f.read())

        #=======================================#
#=======#                NETWORK                #=======#
        #=======================================#

TABLE_REGEX = r'<td [^>]+>\s*%d<br\s*/>\s*<a href="/data/files/([^">]+)[^>]*">'

def save_mos_ecp(n, path):
	logging.debug("Downloading GRECP for element %d" % n)
	logging.debug("Loading table form qchem.pnpi.spb.ru")
	filename = None
	try:
		with urlopen("http://qchem.pnpi.spb.ru/recp") as conn:
			html = conn.read().decode('utf-8')
			m = re.search(TABLE_REGEX % n, html, re.DOTALL)
			halt(not m, "...ELEMENT %d NOT FOUND IN TABLE" % n)
			filename = m.group(1)
			url = 'http://qchem.pnpi.spb.ru/data/files/' + filename
			logging.debug("GRECP url found: %s" % url)
			with urlopen(url) as conn:
				data = conn.read()
				if filename.endswith('.inp'):
					logging.debug("Direct inp link.")
					with open(path, 'w') as ecpf:
						ecpf.write(data.decode('utf-8'))
				elif filename.endswith('.zip'):
					logging.debug("Zip archive. Looking into...")
					file_like_object = io.BytesIO(data)
					with  zipfile.ZipFile(file_like_object) as zf:
						candidates = []
						for name in zf.namelist():
							if m := re.match(r'^%s\d+.inp$' % ELEMENTS[n], name):
								candidates.append(name)
						logging.debug("Inps found: %s" % candidates)
						if not candidates:
							halt(True, "NO GRECP FILES IN ARCHIVE")
						candidate = None
						if len(candidates) == 1:
							candidate = candidates[0]
						else:
							_, candidate = multiselect(' Multiple GRECP files found. Please make a choice:', candidates)
						if candidate:
							with open(path, 'w') as ecpf:
								with zf.open(candidate) as inpf:
									ecpf.write(inpf.read().decode('utf-8'))
				else:
					halt(True, "UNRECOGNIZED FORMAT: %s" % filename)
	except Exception as e:
		halt(True, "DOWNLOAD ERROR: %s" % str(e))

def get_bse_basises(n):
	try:
		url = 'https://www.basissetexchange.org/api/metadata/'
		with urlopen(url) as conn:
			data = conn.read().decode('utf-8')
			j = json.loads(data)
			names = []
			for k in j.keys():
				if 'gto' in j[k]['function_types'] or 'gto_spherical' in j[k]['function_types']:
					v = sorted(j[k]['versions'])[-1]
					if str(n) in j[k]['versions'][v]['elements']:
						names.append(k)
			for name in names:
				logging.info(name)
	except Exception as e:
		halt(True, "ERROR WHILE GETTING BASIS LIST: %s" % str(e))

def get_bse_basis(n, basis_name):
	url = 'https://www.basissetexchange.org/api/basis/%s/format/json/?elements=%d' % (quote(basis_name), n)
	logging.debug(url)
	try:
		with urlopen(url) as conn:
			data = conn.read().decode('utf-8')
			j = json.loads(data)
			if 'error' in j.keys():
				halt(True, "SELECTED BASIS UNAVAILABLE FOR SELECTED ELEMENT")
			if 'elements' not in j.keys() or str(n) not in j['elements'].keys():
				halt(True, "UNEXPECTED ERROR")
			if 'electron_shells' not in j['elements'][str(n)].keys():
				halt(True, "NO FUNCTIONS IN SELECTED BASIS")
			functions = {}
			for shell in j['elements'][str(n)]['electron_shells']:
				for l in shell['angular_momentum']:
					if l not in functions.keys():
						functions[l] = []
					for coef_array in shell['coefficients']:
						cgf = GaussFunctionContracted([])
						for e, c in zip(shell['exponents'], coef_array):
							gf = GaussFunctionNormed(float(e), l)
							cgf.primitives.append([float(c), gf])
						functions[l].append(cgf)
			return Basis(n, functions)
	except Exception as e:
		halt(True, "ERROR WHILE GETTING BASIS: %s" % str(e))

def save_basis_callback(n, name, path):
	try:
		b = get_bse_basis(n, name)
		if b:
			write_basis_nw(b, path)
			return True
	except Exception as e:
		return False

def get_bse_basises_callback(n, iv):
	try:
		get_bse_basises(n)
		iv.active = False
		return False
	except Exception as e:
		return False

def save_starting_basis(n, path):
	default_basis = None
	if n < 87:
		default_basis = 'def2-tzvp'
	elif n > 89 and n < 104:
		default_basis = 'def2-mtzvp'
	variants = []
	if default_basis:
		variants.append(InputVariant("Press ENTER to try %s" % default_basis, r"^$", lambda s, _: save_basis_callback(n, default_basis, path)))
	variants.append(InputVariant("Enter other basis name to try it", r"^[\w+\-(),. ][\w!+\-(),. ]*$", lambda s, _: save_basis_callback(n, s, path)))
	variants.append(InputVariant("Enter '!list' to list all basis sets", r"^\!list$", lambda s, iv: get_bse_basises_callback(n, iv)))
	complexinput('Select basis to get exponents from.', variants)

def steal_polarization(ecp, basis):
	single_basis = None
	double_basis = None
	if ecp.n < 87:
		single_basis = 'def2-tzvp'
		double_basis = 'def2-tzvpp'
	elif ecp.n == 89:
		double_basis = 'stuttgart rlc'
	elif ecp.n > 89 and ecp.n < 104:
		single_basis = 'def2-mtzvp'
		double_basis = 'stuttgart rlc'
	variants = []
	basises = []
	if single_basis:
		variants.append('Steal from %s (one or more functions)' % single_basis)
		basises.append(single_basis)
	if double_basis:
		variants.append('Steal from %s (two or more functions)' % double_basis)
		basises.append(double_basis)
	if variants:
		variants.append('Thou shalt not steal!')
		basises.append(None)
		desired_maxl = estimate_maxl(ecp.n, ecp.ncore)
		desired_maxl = max(desired_maxl, basis.maxl())
		index, variant = multiselect("Steal polarization functions for l>%s?" % desired_maxl, variants)
		basis_to_steal = basises[index]
		if basis_to_steal:
			try:
				b = get_bse_basis(ecp.n, basis_to_steal)
				for l in b.functions.keys():
					if l > desired_maxl:
						basis.functions[l] = b.functions[l]
			except Exception as e:
				logging.warning('The theft has failed.')


        #=======================================#
#=======#          PYSCF CALCULATIONS           #=======#
        #=======================================#


INDEPEND_CALC = 0
OPT_MAIN_CALC = 1
OPT_GRAD_CALC = 2

def check_path(dirname):
	return os.path.join(dirname, 'test.chk')

def calculate_atom(ecp, basis, state, dirname, mode):
	halt(ecp.n != basis.n, "DIFFERENT ATOMIC NUMBERS FOR ECP AND BASIS: %d and %d" % (ecp.n, basis.n))
	element = ELEMENTS[ecp.n]
	mol = pyscf.gto.M(
		verbose = 0,
		atom = '%s 0 0 0' % element,
		basis = {element: basis_to_nw(basis)},
		ecp = {element: ecp_to_pyscf(ecp)},
		charge = state.charge,
		spin = state.spin
	)
	mf = None
	cpath = check_path(dirname)
	if mode == OPT_MAIN_CALC or mode == OPT_MAIN_CALC:
		if os.path.exists(cpath):
			mf = pyscf.scf.UHF(mol)
			chk = pyscf.lib.chkfile.load(cpath, 'scf')
			mf.__dict__.update(chk)
			if mode == OPT_MAIN_CALC:
				mf.chkfile = cpath
			mf.run()
		else:
			mf = pyscf.scf.ROHF(mol)
			mf.run()
			dm1 = mf.make_rdm1()
			mf = pyscf.scf.UHF(mol)
			if mode == OPT_MAIN_CALC:
				mf.chkfile = cpath
			mf.run(dm0=dm1)
	else:
		mf = pyscf.scf.ROHF(mol)
		mf.run()
		dm1 = mf.make_rdm1()
		mf = pyscf.scf.UHF(mol)
		mf.run(dm0=dm1)
	return mf

class AtomicOrbital:
	def __init__(self, coeffs, energy, occ):
		self.coeffs = coeffs
		self.energy = energy
		self.occ = occ

	def l(self):
		lmax = 0
		cmax = 0
		for l in self.coeffs.keys():
			ccur = 0
			for n in self.coeffs[l].keys():
				for c in self.coeffs[l][n]:
					ccur += c ** 2
			if ccur > cmax:
				cmax = ccur
				lmax = l
		return lmax


def read_orbitals(exps, res):
	orbitals = []
	for n in range(res.mo_coeff.shape[-1]):
		for s in range(2):
			coeffs = {}
			energy = res.mo_energy[s][n]
			occ = res.mo_occ[s][n]
			row = 0
			for l in sorted(exps.functions.keys()):
				coeffs[l] = {}
				for i, f in enumerate(exps.functions[l]):
					exp_coeffs = []
					for m in range(2 * l + 1):
						exp_coeffs.append(res.mo_coeff[s][row + m][n])
					row += 2 * l + 1
					coeffs[l][i] = exp_coeffs
			orb = AtomicOrbital(coeffs, energy, occ)
			orbitals.append(orb)
	orbitals = sorted(orbitals, key=lambda o: o.energy)
	for o in orbitals:
		if o.occ:
			logging.debug("%s %8.3f" % (SHELLS[o.l()], o.energy))
	return orbitals

def build_shells(orbitals):
	atom_shells = {}
	for o in orbitals:
		if o.occ > 0:
			l = o.l()
			if not l in atom_shells:
				atom_shells[l] = [[o]]
			else:
				if len(atom_shells[l][-1]) == 4 * l + 2:
					atom_shells[l].append([o])
				else:
					atom_shells[l][-1].append(o)
	return atom_shells


        #=======================================#
#=======#             EXPONENT SET              #=======#
        #=======================================#

class Contraction:
	def __init__(self, l, startn, coeffs):
		self.l = l
		self.startn = startn
		self.coeffs = coeffs

	def nset(self):
		return set(range(self.startn, self.startn + len(self.coeffs)))


class ExponentSet:
	def __init__(self, n, functions):
		self.n = n
		self.functions = functions #{l: [e, opt]}
		self.contractions = {}

	def add_contraction(self, c):
		if c.l not in self.contractions.keys():
			self.contractions[c.l] = [c]
		else:
			for c1 in self.contractions[c.l]:
				halt(c.nset() & c1.nset(), "OVERLLAPPING CONTRACTIONS FOUND!!!")
			self.contractions[c.l].append(c)

	def build_basis(self):
		new_functions = {}
		for l in self.functions.keys():
			new_functions[l] = []
			for n, (e,_) in enumerate(self.functions[l]):
				contracted = False
				if l in self.contractions.keys():
					for c in self.contractions[l]:
						if c.startn == n:
							primitives = []
							for i, coeff in enumerate(c.coeffs):
								primitives.append([coeff, GaussFunctionNormed(self.functions[l][n + i][0], l)])
							cgf = GaussFunctionContracted(primitives)
							cgf.normalize()
							new_functions[l].append(cgf)
							contracted = True
							break
						elif n in c.nset():
							contracted = True
							break
				if not contracted:
					cgf = GaussFunctionContracted([[1, GaussFunctionNormed(e, l)]])
					new_functions[l].append(cgf)
		return Basis(self.n, new_functions)

	def write_to_file(self, filename):
		with open(filename, 'w') as f:
			f.write("%s\n" % ELEMENTS[self.n])
			for l in sorted(self.functions.keys()):
				f.write("%s\n" % SHELLS[l])
				for e, opt in self.functions[l]:
					line = "%20.8f" % e
					if opt:
						line += ' *'
					line += '\n'
					f.write(line)

	def maxl(self):
		return max(self.functions.keys())

	def as_list(self):
		res = []
		for l in sorted(self.functions.keys()):
			for e, o in self.functions[l]:
				res.append((l, e, o))
		return res

	@staticmethod
	def from_list(n, ls):
		functions = {}
		for l, e, o in ls:
			if not (l in functions.keys()):
				functions[l] = []
			functions[l].append([e,o])
		return ExponentSet(n, functions)

	def get_variables(self, l):
		if not (l in self.functions.keys()):
			return []
		xs = []
		for e, o in self.functions[l]:
			if o:
				xs.append(math.log(e))
		return xs

	def set_variables(self, l, xs):
		halt(not (l in self.functions.keys()), "WRONG L")
		for eo in self.functions[l]:
			e, o = eo
			if o:
				halt(len(xs) == 0, "NOT ENOUGH PARAMETERS TO SET")
				eo[0] = math.exp(xs[0])
				xs = xs[1:]
		halt(len(xs) > 0, "TOO MANY PARAMETERS TO SET")

	def get_constraints(self, l, settings):
		if not (l in self.functions.keys()):
			return []
		cons = []
		d1, d2 = settings.get_constraints(l)
		pos = 0
		for i in range(len(self.functions[l]) - 1):
			e1, o1 = self.functions[l][i]
			e2, o2 = self.functions[l][i + 1]
			if o1 and o2:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d1=d1: x[pos] - x[pos+1] - math.log(d1)})
				pos += 1
			elif o1:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d1=d1,e2=e2: x[pos] - math.log(e2) - math.log(d1)})
				pos += 1
			elif o2:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d1=d1,e1=e1: math.log(e1) - x[pos+1] - math.log(d1)})
		pos = 0
		for i in range(len(self.functions[l]) - 2):
			e1, o1 = self.functions[l][i]
			e2, o2 = self.functions[l][i + 1]
			e3, o3 = self.functions[l][i + 2]
			if o1 and o2 and o3:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d2=d2: x[pos] - x[pos+2] - math.log(d2)})
				pos += 1
			elif o1 and o3:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d2=d2: x[pos] - x[pos+1] - math.log(d2)})
				pos += 1
			elif o1:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d2=d2,e3=e3: x[pos] - math.log(e3) - math.log(d2)})
				pos += 1
			elif o3:
				cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos,d2=d2,e1=e1: math.log(e1) - x[pos+1] - math.log(d2)})
		return cons

	@staticmethod
	def from_basis(basis):
		functions = {}
		for l in basis.functions.keys():
			exponents = set()
			for cgf in basis.functions[l]:
				for _, gf in cgf.primitives:
					exponents.add(gf.a)
			functions[l] = [[x, False] for x in sorted(list(exponents), reverse=True)]
		return ExponentSet(basis.n, functions)

	@staticmethod
	def from_file(filename):
		functions = {}
		with open(filename) as f:
			lines = f.read().splitlines()
			_, n = get_element(lines[0].strip())
			l = -1
			for i, line in enumerate(lines[1:]):
				if "contraction" in line.lower():
					i -= 1
					break
				if m := re.match(r'^[spdfghi]$', line.strip()):
					l = SHELLS.index(m.group())
					functions[l] = []
					continue
				halt(l == -1, 'EXPONENTS INPUT ERROR')
				try:
					ls = line.split()
					e = float(ls[0])
					opt = False
					if len(ls) > 1:
						halt((len(ls) > 2) or (ls[1] != '*'), "CANNOT PARSE EXPONENTS LINE: '%s'" % line)
						opt = True
					functions[l].append([e, opt])
				except Exception:
					halt(True, "CANNOT PARSE EXPONENTS LINE: '%s'" % line)
			return ExponentSet(n, functions)


def set_radii(exps, radii):
	for l in sorted(exps.functions.keys()):
		if l in radii.keys():
			r = max(radii[l], 1e-5)
			a = (l + 1) / r**2 / 2
			logging.debug("L=%d; r=%.2f; a=%.2f" % (l, r, a))
			for eo in exps.functions[l]:
				if eo[0] >= a:
					eo[1] = True

def set_single_radius_callback(exps, s, desired_maxl):
	r = float(s)
	radii = {l: r for l in range(desired_maxl + 1)}
	set_radii(exps, radii)
	return True

def set_multi_radii_callback(exps, s, desired_maxl):
	rs = [float(x) for x in re.split(r'[, ]+', s)]
	if len(rs) != desired_maxl + 1:
		logging.warning('Wrong number of radii')
		return False
	radii = {l: r for l, r in zip(range(desired_maxl + 1), rs)}
	set_radii(exps, radii)
	return True

def print_vdz_radii_callback(exps, n, ncore, iv):
	if iv.custom_field:
		set_radii(exps, iv.custom_field)
		return True
	try:
		arv_basis = get_bse_basis(n, "ano-rcc-vdz")
		radii = estimate_core_bounds(arv_basis, ncore)
		print_subtitle('Estimated radii')
		for l in radii.keys():
			r = radii[l]
			a = (l + 1) / r**2 / 2
			logging.info("L=%d; r=%.2f; a=%.2f" % (l, r, a))
			logging.info('')
		iv.custom_field = radii
		iv.regex = "^$"
		iv.hint = "Press ENTER to use estimated ANO-RCC-VDZ radii"
		return False
	except Exception:
		return False


def select_vc_divide(exps, n, ncore, desired_maxl):
	float_regex = r'(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?'
	bounds = None
	print_subtitle("Valence/core divide selection")
	variants = []
	variants.append(InputVariant("Enter single radius for all Ls", "^%s$" % float_regex, lambda s, _: set_single_radius_callback(exps, s, desired_maxl)))
	variants.append(InputVariant("Enter radii sequence for Ls in 0-%d range" % desired_maxl, r"^%s$" % r'[, ]+'.join([float_regex] * (desired_maxl + 1)), lambda s, _: set_multi_radii_callback(exps, s, desired_maxl)))
	variants.append(InputVariant("Enter '!vdz' to print max radii of ANO-RCC-VDZ functions", r"^\!vdz$", lambda s, iv: print_vdz_radii_callback(exps, n, ncore, iv)))
	complexinput('Select how to divide exponents:', variants)


        #=======================================#
#=======#        EXPONENT OPTIMIZATION          #=======#
        #=======================================#


class ExpOptData:
	def __init__(self, exps):
		self.start_exps = ExponentSet.from_list(exps.n, exps.as_list())
		self.current_exps = ExponentSet.from_list(exps.n, exps.as_list())
		self.estart = None
		self.emin = None
		self.elast = None

	def copy_start_exps(self):
		return ExponentSet.from_list(self.start_exps.n, self.start_exps.as_list())

	def update_from(self, exps, l):
		xs = exps.get_variables(l)
		self.current_exps.set_variables(l, xs)

	def validate_result(self, ecur):
		if (self.elast is None) or (self.emin is None):
			return True
		return ecur <= (self.elast + self.emin) / 2 + 1e-5

	def process_result(self, ecur, restart_path):
		self.elast = ecur
		if (self.estart is None):
			self.estart = ecur
		if (self.emin is None) or (ecur <= self.emin):
			self.current_exps.write_to_file(restart_path)
			self.emin = ecur


def calculate_exp(xs, l, exps, ecp, state, eps, restart_path, dirname, expdata):
	exps.set_variables(l, xs)
	expdata.update_from(exps, l)
	t = check_if_step_not_grad(xs, eps)
	if t:
		logging.debug('%s' % [math.exp(x) for x in xs])
	basis = exps.build_basis()
	res = calculate_atom(ecp, basis, state, dirname, OPT_MAIN_CALC if t else OPT_GRAD_CALC).e_tot
	if t:
		logging.debug('%f' % res)
		if not expdata.validate_result(res):
			logging.debug("Going up: resetting checkfile")
			cpath = check_path(dirname)
			if os.path.exists(cpath):
				os.remove(cpath)
			res = calculate_atom(ecp, basis, state, dirname, OPT_MAIN_CALC).e_tot
			logging.debug('%f' % res)
		expdata.process_result(res, restart_path)
	return res


def optimize_exponents(exps, ecp, state, cs, maxit, dirname, restart_filename, lset):
	opt_delta_total = 0
	restart_path = os.path.join(dirname, restart_filename)
	cpath = check_path(dirname)
	if os.path.exists(cpath):
		os.remove(cpath)
	Estart = calculate_atom(ecp, exps.build_basis(), state, dirname, INDEPEND_CALC).e_tot
	logging.info("Non-optimized E=%.6f a.u." % (Estart))
	expdata = ExpOptData(exps)
	lrange = sorted(lset) if lset else range(exps.maxl() + 1)
	for l in lrange:
		exps_optimizing = expdata.copy_start_exps()
		expdata.estart = None
		expdata.emin = None
		opt_delta_l = 0
		print_subtitle("Optimizing %s-component" % SHELLS[l].upper())
		x0 = exps_optimizing.get_variables(l)
		if not x0:
			continue
		cons = exps_optimizing.get_constraints(l, cs)
		optsteps = generate_optsteps(cs.get_constraints(l)[0])
		global last_xe_step
		global last_xe
		last_xe_step = None
		last_xe = None
		for i, (eps, ftol) in enumerate(optsteps):
			expdata.estart = None
			x0 = exps_optimizing.get_variables(l)
			logging.info("%s-optimization step %d: eps=%.2e, ftol=%.2e" % (SHELLS[l].upper(), i+1, eps, ftol))
			f = lambda x: calculate_exp(x, l, exps_optimizing, ecp, state, eps, restart_path, dirname, expdata)
			max_bound = max(x0) + math.log(10)
			min_bound = min(math.log(1), min(x0) - math.log(2))
			bounds = [(min_bound,max_bound)] * len(x0)
			res = scipy.optimize.minimize(f, x0, args=(), method='SLSQP', jac=None, bounds=bounds, constraints=cons, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
			if res.success:
				logging.info("... finished in %d iterations and %d evaluations." % (res.nit, res.nfev))
				de_step = res.fun - expdata.estart
				logging.info("... %s-step %d: dE=%.6f a.u." % (SHELLS[l].upper(), i+1, de_step))
				opt_delta_l += de_step
			else:
				halt(True, 'OPTIMIZATION FAILED WITH STATUS %d' % res.status)
			exps_optimizing.set_variables(l, res.x)
			expdata.update_from(exps_optimizing, l)
		opt_delta_total += opt_delta_l
		logging.info("Optimization for %s-component finished" % SHELLS[l])
		logging.info("%s-component dE=%.6f a.u." % (SHELLS[l], opt_delta_l))
	logging.info("Basis optimization finished")
	logging.info("Sum of independent dEs=%.6f a.u." % opt_delta_total)
	exps = expdata.current_exps
	Efinish = calculate_atom(ecp, exps.build_basis(), state, dirname, INDEPEND_CALC).e_tot
	exps.write_to_file(restart_path)
	logging.info("Optimized E=%.6f a.u." % (Efinish))
	logging.info("Total dEs=%.6f a.u." % (Efinish - Estart))
	return exps

class ConstraintSettings:

	def __init__(self, maxl, constraints):
		self.maxl = maxl
		self.constraints = constraints #{l: [c1, c2]}

	def default_constraints(self, l):
		if l == 0:
			return [1.3, 1.7]
		elif l == 1:
			return [1.2, 1.5]
		else:
			return [1.1, 1.3]

	def get_constraints(self, l):
		if l in self.constraints.keys():
			return self.constraints[l]
		return self.default_constraints(l)

	def write_to_file(self, filename):
		with open(filename, 'w') as f:
			for l in range(self.maxl + 1):
				f.write('%s: %4.1f, %4.1f\n' % tuple([SHELLS[l]] + self.get_constraints(l)))

	@staticmethod
	def from_file(filename):
		constraints = {}
		with open(filename) as f:
			for line in f.read().splitlines():
				if m := re.match(r'^([spdfghi]):\s+([\d.]+),\s*([\d.]+)+$', line):
					l = SHELLS.index(m.group(1))
					c1 = float(m.group(2))
					c2 = float(m.group(3))
					constraints[l] = [c1, c2]
				else:
					halt(True, "CANNOT PARSE CONSTRAINTS LINE: '%s'" % line)
		return ConstraintSettings(max(constraints.keys()), constraints)

def generate_optsteps(constraint):
	tol = 1e-4
	optsteps = []
	for i in range(3):
		optsteps += [[math.log(constraint) * 0.5 / 10**i, tol / 10**i]]
	return optsteps


        #=======================================#
#=======#       CONTRACTION OPTIMIZATION        #=======#
        #=======================================#

CDAMP = 3
EDAMP = 1

def get_damped_variables(cgf):
	xc, index = cgf.get_coeff_variables()
	xc = [x * CDAMP for x in xc]
	xe = cgf.get_exp_variables()
	xe = [x * EDAMP for x in xe]
	return xc + xe, index

def get_damped_bounds(cgf):
	xc, index = cgf.get_coeff_variables()
	xc = [x * CDAMP for x in xc]
	xe = cgf.get_exp_variables()
	xe = [x * EDAMP for x in xe]
	bc = [(-10 * CDAMP, 10 * CDAMP)] * len(xc)
	be = [(min(xe) - EDAMP, max(xe) + EDAMP)] * len(xe)
	return bc + be

def get_damped_constraints(cgf):
	xc, index = cgf.get_coeff_variables()
	xc = [x * CDAMP for x in xc]
	xe = cgf.get_exp_variables()
	xe = [x * EDAMP for x in xe]
	pos = len(xc)
	cons = []
	for i in range(len(xe) - 1):
		cons.append({'type': 'ineq', 'fun' : lambda x,pos=pos: x[pos] - x[pos+1]})
		pos += 1
	return cons

def set_damped_variables(cgf, xd, index):
	halt(len(xd) != (len(cgf.primitives) * 2 - 1), "SOMETHING VERY WRONG HAPPENED!")
	xc = [x / CDAMP for x in xd[:len(xd) // 2]]
	xe = [x / EDAMP for x in xd[len(xd) // 2:]]
	cgf.set_coeff_variables(xc, index)
	cgf.set_exp_variables(xe)

def log_damped_variables(xd):
	xc = [x / CDAMP for x in xd[:len(xd) // 2]]
	xe = [math.exp(x / EDAMP) for x in xd[len(xd) // 2:]]
	logging.debug('%s' % (list(xc)))
	logging.debug('%s' % (list(xe)))

class ContrOptData:
	def __init__(self, cgf):
		self.cgf = cgf
		self.estart = None
		self.emin = None
		self.elast = None

	def validate_result(self, ecur):
		if (self.elast is None) or (self.emin is None):
			return True
		return ecur <= (self.elast + self.emin) / 2 + 1e-5

	def process_result(self, ecur, basis, restart_path):
		self.elast = ecur
		if (self.estart is None):
			self.estart = ecur
		if (self.emin is None) or (ecur <= self.emin):
			write_basis_nw(basis, restart_path)
			self.emin = ecur

def calculate_contr(xs, basis, coptdata, index, ecp, state, eps, restart_path, dirname):
	set_damped_variables(coptdata.cgf, xs, index)
	t = check_if_step_not_grad(xs, eps)
	if t:
		log_damped_variables(xs)
	res = calculate_atom(ecp, basis, state, dirname, OPT_MAIN_CALC if t else OPT_GRAD_CALC).e_tot
	if t:
		logging.debug('%f' % res)
		if not coptdata.validate_result(res):
			logging.debug("Going up: resetting checkfile")
			cpath = check_path(dirname)
			if os.path.exists(cpath):
				os.remove(cpath)
			res = calculate_atom(ecp, basis, state, dirname, OPT_MAIN_CALC).e_tot
			logging.debug('%f' % res)
		coptdata.process_result(res, basis, restart_path)
	return res

def optimize_contraction(basis, ecp, state, maxit, dirname, restart_filename):
	opt_delta_total = 0
	restart_path = os.path.join(dirname, restart_filename)
	cpath = check_path(dirname)
	if os.path.exists(cpath):
		os.remove(cpath)
	Estart = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC).e_tot
	logging.info("Non-optimized E=%.6f a.u." % (Estart))
	for l in sorted(basis.functions.keys()):
		for cgf in basis.functions[l]:
			if len(cgf.primitives) > 1:
				print_subtitle("Current contraction to optimize (L=%d):" % l)
				for c, f in cgf.primitives:
					logging.info("%12.6f %12.6f" % (f.a, c))
				optimize = yesno("Optimize?")
				if optimize:
					optsteps = [(0.001, 1e-4), (0.0003, 1e-5), (0.0001, 1e-6)]
					global last_xe_step
					global last_xe
					last_xe_step = None
					last_xe = None
					opt_delta_f = 0
					for i, (eps, ftol) in enumerate(optsteps):
						x0, index = get_damped_variables(cgf)
						coptdata = ContrOptData(cgf)
						logging.info("Optimization step %d: eps=%.2e, ftol=%.2e" % (i+1, eps, ftol))
						bounds = get_damped_bounds(cgf)
						constraints = get_damped_constraints(cgf)
						f = lambda x: calculate_contr(x, basis, coptdata, index, ecp, state, eps, restart_path, dirname)
						res = scipy.optimize.minimize(f, x0, args=(), method='SLSQP', jac=None, bounds=bounds, constraints=constraints, tol=None, callback=None, options={'disp': False, 'eps': eps, 'maxiter': maxit, 'ftol': ftol})
						if res.success:
							logging.info("... finished in %d iterations and %d evaluations." % (res.nit, res.nfev))
							logging.info("... step %d: dE=%.6f a.u." % (i+1, res.fun - coptdata.estart))
							opt_delta_f += res.fun - coptdata.estart
						else:
							halt(True, 'OPTIMIZATION FAILED WITH STATUS %d' % res.status)
						set_damped_variables(cgf, res.x, index)
					opt_delta_total += opt_delta_f
					logging.info("Optimization finished: dE=%.6f a.u." % (opt_delta_f))
	logging.info("Contraction optimization finished")
	logging.info("Contraction dE=%.6f a.u." % opt_delta_total)
	Efinish = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC).e_tot
	logging.info("Optimized E=%.6f a.u." % (Efinish))
	logging.info("Total dEs=%.6f a.u." % (Efinish - Estart))
	return basis


        #=======================================#
#=======#              ATOM STATE               #=======#
        #=======================================#

class AtomState:
	def __init__(self, charge, spin, comment):
		self.charge = charge
		self.spin = spin
		self.comment = comment

	def write_to_file(self, filename):
		with open(filename, 'w') as f:
			f.write('charge = %d\n' % self.charge)
			f.write('spin = %d\n' % self.spin)
			f.write('comment = %s\n' % self.comment)

	@staticmethod
	def from_file(filename):
		with open(filename) as f:
			lines = f.read().splitlines()
			m = re.match(r"^charge\s*=\s*(-?\d+)$", lines[0].strip())
			halt(not m, "CANNOT PARSE ATOMSTATE LINE: %s" % lines[0])
			charge = int(m.group(1))
			m = re.match(r"^spin\s*=\s*(\d+)$", lines[1].strip())
			halt(not m, "CANNOT PARSE ATOMSTATE LINE: %s" % lines[1])
			spin = int(m.group(1))
			m = re.match(r"^comment\s*=\s*(.+)$", lines[2].strip())
			halt(not m, "CANNOT PARSE ATOMSTATE LINE: %s" % lines[2])
			comment = m.group(1)
			return AtomState(charge, spin, comment)

def unquote(nist_csv_cell):
	regex = r'''^="(.+?)"$'''
	if m := re.match(regex, nist_csv_cell):
		return m.group(1)
	return nist_csv_cell

def get_nist_states(n, val):
	logging.debug("Downloading ion shells for element %d from physics.nist.gov" % n)
	url = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=%s&units=1&format=2&remove_js=on&order=0&ion_charge_out=on&shells_out=on&level_out=on&e_out=0&submit=Retrieve+Data' % (ELEMENTS[n].lower())
	logging.debug(url)
	try:
		with urlopen(url) as conn:
			data = conn.read().decode('utf-8')
			states = []
			csv_input = csv.reader(io.StringIO(data), delimiter=',')
			rows = list(csv_input)
			headers = rows[0]
			charge_index = headers.index('Ion Charge')
			conf_index = headers.index('Ground Shells (a)') if 'Ground Shells (a)' in headers else headers.index('Ground Shells')
			term_index = headers.index('Ground Level')
			core = None
			for row in rows[1:]:
				if len(row) != len(headers):
					break
				row = [unquote(x) for x in row]
				charge = int(row[charge_index])
				if charge > val:
					break
				m = re.match(r'^\d+', row[term_index])
				if not m:
					continue
				spin = int(m.group()) - 1
				conf = row[conf_index]
				m = re.match(r'^\[\w+\]', conf)
				conf_core = None
				if m:
					conf_core = m.group()
				else:
					conf_core = ''
				if core is None:
					core = conf_core
				if core != conf_core:
					break
				state = AtomState(charge, spin, 'NIST data: %s' % conf)
				states.append(state)
	except Exception as e:
		halt(True, "ERROR WHILE GETTING CONFIGURATIONS: %s" % str(e))
	return states

def get_nist_configuration(n):
	logging.debug("Downloading configuration for element %d from physics.nist.gov" % n)
	url = 'https://physics.nist.gov/cgi-bin/ASD/ie.pl?spectra=%s&units=1&format=2&remove_js=on&order=0&ion_charge_out=on&shells_out=on&level_out=on&e_out=0&submit=Retrieve+Data' % (ELEMENTS[n].lower())
	logging.debug(url)
	try:
		with urlopen(url) as conn:
			data = conn.read().decode('utf-8')
			states = []
			csv_input = csv.reader(io.StringIO(data), delimiter=',')
			rows = list(csv_input)
			headers = rows[0]
			conf_index = headers.index('Ground Shells (a)')
			row = rows[1]
			if len(row) != len(headers):
				return None
			row = [unquote(x) for x in row]
			conf = row[conf_index]
			m = re.match(r'^(?:\[\w+\]\.)?(\d[spdfg]\d*(\.\d[spdfg]\d*)*)$', conf)
			configuration = {}#n,l : occ
			s = m.group(1)
			for part in s.split('.'):
				m = re.match(r'^(\d)([spdfg])(\d*)$', part)
				n = int(m.group(1))
				l = int(SHELLS.index(m.group(2)))
				occ = int(m.group(3)) if m.group(3) else 1
				configuration[(n,l)] = occ
			return configuration
	except Exception as e:
		halt(True, "ERROR WHILE GETTING CONFIGURATIONS: %s" % str(e))

def select_charge_callback(s, states, filepath):
	i = int(s)
	if i < len(states):
		state = states[i]
		state.write_to_file(filepath)
		logging.debug('Selected state: c=%d, s=%d, %s' % (state.charge, state.spin, state.comment))
		return True
	else:
		logging.warning('Charge %d is not in list.' % i)
		return False

def select_charge_spin_callback(s, val, filepath):
	m = re.match(r'^(-\d+),?\s*(\d+)$', s)
	c = int(m.group(1))
	s = int(m.group(2))
	if c < val + 1:
		#TODO: validate spin value
		state = AtomState(c, s, 'User-defined')
		state.write_to_file(filepath)
		logging.debug('Selected state: c=%d, s=%d, %s' % (state.charge, state.spin, state.comment))
		return True
	else:
		logging.warning('Charge %d is too large.' % c)
		return False

def get_nist_states_callback(n, states, val, iv):
	try:
		nist_states = get_nist_states(n, val)
		states.clear()
		states.extend(nist_states)
		print_subtitle("NIST atom states:")
		for state in states:
			logging.info(" charge=%d; spin=%d (%s)" % (state.charge, state.spin, state.comment))
		iv.active = False
		return False
	except Exception as e:
		raise e
		return False

def select_atom_state(n, nval, filepath):
	val, conf_s, spin = estimate_ion_valence_shells(n)
	val = min(val, nval)
	states = [AtomState(0, spin, 'Auto-estimated as %s' % conf_s)]
	for i in range(val):
		_, conf_s1, spin1 = estimate_ion_valence_shells(n, i+1)
		states.append(AtomState(i + 1, spin1, 'Auto-estimated as %s' % conf_s1))
	print_subtitle("Estimated atom states:")
	for state in states:
		logging.info(" charge=%d; spin=%d (%s)" % (state.charge, state.spin, state.comment))
	variants = []
	variants.append(InputVariant("Enter charge to select state from list", r"^\d+$", lambda s, _: select_charge_callback(s, states, filepath)))
	variants.append(InputVariant("Enter charge, spin to define the state", r"^(-\d+),?\s*(\d+)$", lambda s, _: select_charge_spin_callback(s, val, filepath)))
	variants.append(InputVariant("Enter '!nist' to update state configuration from NIST", r"^\!nist$", lambda _, iv: get_nist_states_callback(n, states, val, iv)))
	complexinput('Select predefined state or input charge and spin.', variants)



        #=======================================#
#=======#              DISCARDING               #=======#
        #=======================================#



def discarded_copy(exps, L, E):
	ls = []
	for l, e, o in exps.as_list():
		if l != L or e != E:
			ls.append((l, e, o))
	return ExponentSet.from_list(exps.n, ls)

def check_for_discarding(exps, ecp, state, dirname):
	LIMIT_MAX = 0.05
	LIMIT_MEAN = 0.005
	print_subtitle('Checking orbital coefficients.')
	basis = exps.build_basis()
	res = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC)
	orbitals = read_orbitals(exps, res)
	exps_to_discard = []
	discarded_ls = set()

	for l in sorted(exps.functions.keys()):
		logging.debug('---------')
		for i, (e, o) in enumerate(exps.functions[l]):
			if o:
				c = []
				for orb in orbitals:
					if orb.occ > 0 and orb.l() == l:
						c.extend([abs(x) for x in orb.coeffs[l][i]])
				if c:
					mx = max(c)
					mn = sum(c) / len(c)
					if mx < LIMIT_MAX and mn < LIMIT_MEAN:
						exps_to_discard.append((l, e, mx, mn))
					logging.debug("L=%d; a=%15.6f; max=%10.2e; mean=%10.2e" % (l, e, mx, mn))

	if exps_to_discard:
		logging.info('Candidates for exclusion found:')
	for l, e, mx, mn in exps_to_discard:
		logging.info("L=%d; a=%15.6f; max=%10.2e; mean=%10.2e" % (l, e, mx, mn))
	for l, e, mx, mn in exps_to_discard:
		variants = ['Discard', 'Keep', 'Check dE']
		index, variant = multiselect("Discard L=%d; a=%.6f; max=%.2e; mean=%.2e?" % (l, e, mx, mn), variants)
		if index == 0:
			exps = discarded_copy(exps, l, e)
			discarded_ls.add(l)
		elif index == 1:
			continue
		else:
			basis_before = exps.build_basis()
			basis_after = discarded_copy(exps, l, e).build_basis()
			e_before = calculate_atom(ecp, basis_before, state, dirname, INDEPEND_CALC).e_tot
			e_after = calculate_atom(ecp, basis_after, state, dirname, INDEPEND_CALC).e_tot
			variants = ['Discard', 'Keep']
			index, variant = multiselect("E(with)=%.6f; E(without)=%.6f; dE=%.6f. Discard?" % (e_before, e_after, e_after - e_before), variants)
			if index == 0:
				exps = discarded_copy(exps, l, e)
				discarded_ls.add(l)
			elif index == 1:
				continue
	return exps, discarded_ls

        #=======================================#
#=======#             CONTRACTING               #=======#
        #=======================================#

def contracted_copy(exps, L, a1, a2, coeff):
	ls = []
	for l, e, o in exps.as_list():
		ls.append((l, e, o))
	new_exps = ExponentSet.from_list(exps.n, ls)
	n1 = -1
	n2 = -1
	for n, (e, o) in enumerate(exps.functions[L]):
		if e == a1:
			n1 = n
		elif e == a2:
			n2 = n
	halt(n1 == -1 or n2 == -1 or n2 - n1 != 1, "SOMETHING WRONG HAPPENED!")
	added = False
	for l in exps.contractions.keys():
		for c in exps.contractions[l]:
			new_c = Contraction(c.l, c.startn, c.coeffs[:])
			if (l == L) and not added:
				if new_c.startn + len(new_c.coeffs) - 1 == n1:
					new_c.coeffs.append(coeff * new_c.coeffs[-1])
					added = True
				elif new_c.startn == n2:
					new_c.coeffs = [new_c.coeffs[0] / coeff] + new_c.coeffs
					new_c.startn = n1
					added = True
			new_exps.add_contraction(new_c)
	if not added:
		new_c = Contraction(L, n1, [1, coeff])
		new_exps.add_contraction(new_c)
	return new_exps


	def add_exponent(self, n, coeff):
		if n == startn + len(self.coeffs):
			self.coeffs.append(coeff * self.coeffs[-1])
			return True
		elif n > startn + len(self.coeffs):
			return False
		else:
			halt(True, "ATTEMPT TO ADD OVERLAPPING CONSTRAINT!")

def estimate_segmented_contraction(exps, ecp, state, dirname):
	TOL = 1e-5
	print_subtitle('Estimating possible contractions.')
	basis = exps.build_basis()
	res = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC)
	orbitals = read_orbitals(exps, res)
	logging.info('Uncontracted E = %.6f' % res.e_tot)
	candidates = []

	for l in sorted(exps.functions.keys()):
		logging.debug('---------')
		for i in range(len(exps.functions[l]) - 1):
			a1 = exps.functions[l][i][0]
			a2 = exps.functions[l][i + 1][0]
			weighted_sum = 0
			sum_of_weights = 0
			values = []
			weights = []
			for orb in orbitals:
				if orb.occ > 0 and orb.energy < 0 and orb.l() == l:
					for c1, c2 in zip(orb.coeffs[l][i], orb.coeffs[l][i + 1]):
						if (abs(c1) > TOL or abs(c2) > TOL):
							value = c2 / c1
							weight = abs(orb.energy) * max(abs(c1), abs(c2))
							values.append(value)
							weights.append(weight)
							weighted_sum += value * weight
							sum_of_weights += weight
			if sum_of_weights > 0:
				weighted_mean = weighted_sum / sum_of_weights
				weighted_square_sum = 0
				for value, weight in zip(values, weights):
					weighted_square_sum += weight * (value - weighted_mean)**2
				weighted_deviaton = (weighted_square_sum * len(values) / (len(values) - 1) / sum_of_weights)**0.5 if len(values) > 1 else 0
				relative_weighted_deviaton = abs(weighted_deviaton / weighted_mean)
				logging.debug("L=%d; a1=%12.6f, a2=%12.6f; coeff=%8.2e; deviation=%3d%%" % (l, a1, a2, weighted_mean, relative_weighted_deviaton * 100))
				candidates.append((l, a1, a2, weighted_mean, relative_weighted_deviaton))
			else:
				logging.debug("L=%d; a1=%12.6f, a2=%12.6f; coeff=***" % (l, a1, a2))

	if candidates:
		logging.info('Candidates for contraction found:')
	for l, a1, a2, weighted_mean, relative_weighted_deviaton in candidates:
		if relative_weighted_deviaton < 0.5:
			variants = ['Contract', 'Skip', 'Check dE']
			index, variant = multiselect("Contract L=%d; a1=%.6f, a2=%.6f; coeff=%.2e; deviation=%d%%?" % (l, a1, a2, weighted_mean, relative_weighted_deviaton * 100), variants)
			if index == 0:
				exps = contracted_copy(exps, l, a1, a2, weighted_mean)
			elif index == 1:
				pass
			else:
				basis_before = exps.build_basis()
				basis_after = contracted_copy(exps, l, a1, a2, weighted_mean).build_basis()
				e_before = calculate_atom(ecp, basis_before, state, dirname, INDEPEND_CALC).e_tot
				e_after = calculate_atom(ecp, basis_after, state, dirname, INDEPEND_CALC).e_tot
				variants = ['Contract', 'Skip']
				index, variant = multiselect("E(uncontr)=%.6f; E(contr)=%.6f; dE=%.6f. Contract?" % (e_before, e_after, e_after - e_before), variants)
				if index == 0:
					exps = contracted_copy(exps, l, a1, a2, weighted_mean)
				elif index == 1:
					pass
	return exps.build_basis()

def estimate_general_contraction(exps, ecp, state, dirname):
	print_subtitle('Estimating shell-based contractions.')
	basis = exps.build_basis()
	res = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC)
	Estart = res.e_tot
	logging.info('Uncontracted E = %.6f' % Estart)
	orbitals = read_orbitals(exps, res)
	shells = build_shells(orbitals)
	functions = {}
	for l in shells.keys():
		functions[l] = []
		for shell in shells[l]:
			primitives = []
			for i, (e,_) in enumerate(exps.functions[l]):
				weighted_sum = 0
				sum_of_weights = 0
				for orb in shell:
					if orb.occ > 0 and orb.energy < 0 and orb.l() == l:
						for c in orb.coeffs[l][i]:
							weight = abs(orb.energy)
							weighted_sum += c * weight
							sum_of_weights += weight
				weighted_mean = weighted_sum / sum_of_weights
				primitives.append([weighted_mean, GaussFunctionNormed(e, l)])
			cgf = GaussFunctionContracted(primitives)
			cgf.normalize()
			functions[l].append(cgf)
	basis = Basis(ecp.n, functions)
	Efinish = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC).e_tot
	logging.info('Contracted E = %.6f' % Efinish)
	logging.info('Contraction dE = %.6f' % (Efinish - Estart))
	uncontract = yesno('Uncontract some diffuse orbitals?')
	if uncontract:
		for l in functions.keys():
			cgf = functions[l][-1]
			amax = max(cgf.primitives, key=itemgetter(0))[1].a
			for _, f in cgf.primitives:
				if f.a < amax:
					functions[l].append(GaussFunctionContracted([[1, GaussFunctionNormed(f.a, l)]]))
		basis = Basis(ecp.n, functions)
		Efinish = calculate_atom(ecp, basis, state, dirname, INDEPEND_CALC).e_tot
		logging.info('New E = %.6f' % Efinish)
		logging.info('New dE = %.6f' % (Efinish - Estart))
	return basis

        #=======================================#
#=======#      CRYSTAL BASIS CONFIGURATION      #=======#
        #=======================================#

def write_crystal_callback(ecp, basis, conf_string, filename):
	logging.debug(conf_string)
	occs = {}
	for s in re.findall(r'(s[12]|p[1-6]|d(10|\d)|f(1[1-4]|\d))', conf_string):
		m = re.match(r'([spdf])(\d+)', s[0])
		l = SHELLS.index(m.group(1))
		occ = int(m.group(2))
		if l not in occs.keys():
			occs[l] = []
		occs[l].append(occ)
	write_ecp_and_basis_crystal(ecp, basis, occs, filename)
	return True

def write_crystal_nist_callback(ecp, basis, filename):
	configuration = estimate_ecp_configuration(ecp.n, ecp.ncore)
	valence_configuration = estimate_valence_configuration(ecp.n)
	valence_configuration_map = {}
	try:
		nist_configuration = get_nist_configuration(ecp.n)
	except Exception as e:
		logging.warning('Failed to get NIST configuration')
		return False
	configuration_map = {}
	for n, l, occ in configuration:
		configuration_map[(n,l)] = occ
	for n, l, occ in valence_configuration:
		configuration_map[(n,l)] = 0
	for n, l in nist_configuration.keys():
		configuration_map[(n,l)] = nist_configuration[(n,l)]
	conf_string = ''
	for n, l in sorted(configuration_map.keys(), key=itemgetter(1,0)):
		occ = configuration_map[(n,l)]
		if occ:
			conf_string += '%s%d' % (SHELLS[l], occ)
	return write_crystal_callback(ecp, basis, conf_string, filename)

def estimate_and_write_crystal(ecp, basis, filename):
	configuration = estimate_ecp_configuration(ecp.n, ecp.ncore)
	if configuration:
		conf_string = ''
		for n, l, occ in sorted(configuration, key=itemgetter(1,0)):
			if occ:
				conf_string += '%s%d' % (SHELLS[l], occ)
		variants = []
		variants.append(InputVariant("Press ENTER to use auto-estimation: %s" % conf_string, r"^$", lambda s, _: write_crystal_callback(ecp, basis, conf_string, filename)))
		variants.append(InputVariant("Enter '!nist' to get configuration from NIST", r"\!nist$", lambda s, _: write_crystal_nist_callback(ecp, basis, filename)))
		variants.append(InputVariant("Enter configuration in 's2s1p6d4' format", r"(s[12]|p[1-6]|d(10|\d)|f(\d|1[1-4]))+$", lambda s, iv: write_crystal_callback(ecp, basis, s, filename)))
		complexinput('Select neutral configuration for CRYSTAL basis format', variants)
	else:
		logging.warning("Couldn't estimate configuration for ncore=%d" % ecp.ncore)

        #=======================================#
#=======#                 STEPS                 #=======#
        #=======================================#

class StepManager:
	def __init__(self, steps):
		self.steps = steps
		for i, s in enumerate(self.steps):
			s.index = i

	def run(self):
		i = 0
		timestamp = 0
		while True:
			if i == len(self.steps):
				logging.info("Finished sucessfully!")
				return
			s = self.steps[i]
			halt(i != s.index, 'SOMEONE MESSED WITH STEPS ORDER!')
			print_title('Step %d: %s' % (s.index, s.title))
			logging.info(s.description)
			prev = self.steps[i - 1] if i else None
			s.prepare(prev)
			if s.can_be_skipped(timestamp):
				repeat = yesno('Step results are aready present. Overwrite?')
				if not repeat:
					s.finalize()
					i += 1
					timestamp = s.lastmodified()
					continue
			elif s.dir_exists() and s.must_be_rewritten(timestamp):
				logging.warning('Old step results are present. All subsequent steps will be deleted!')
				logging.warning('Backup everything you want and press ENTER to continue.')
				input()
				for s1 in self.steps[i:]:
					shutil.rmtree(s1.dirname(), ignore_errors=True)
			try:
				s.run()
				timestamp = s.lastmodified()
				i += 1
			except Exception as e:
				repeat = yesno('Step exited with error. Repeat the step?')
				if repeat:
					continue
				else:
					raise e
			s.finalize()

class Step:
	def __init__(self):
		self.index = 0
		self.title = 'Empty'
		self.description = 'Doing nothing'
		self.input_files = []
		self.output_files = []
		self.prev_class = type(None)
		self.prev = None
		self.ecp = None
		self.basis = None

	def dirname(self):
		return 'Step %d - %s' % (self.index, self.title)

	def lastmodified(self):
		halt(not self.output_files, 'NO OUTPUT FILES SET FOR TASK. CANNOT GET MODIFIED TIME!')
		halt(not self.output_files_exist(), 'NOT ALL OUTPUT FILES EXIST. CANNOT GET MODIFIED TIME!')
		f = self.output_files[0]
		t = os.path.getmtime(os.path.join(self.dirname(), f))
		for f in self.output_files[1:]:
			t = max(t, os.path.getmtime(os.path.join(self.dirname(), f)))
		return t

	def firstmodified(self):
		files = self.input_files + self.output_files
		files = [f for f in files if os.path.exists(os.path.join(self.dirname(), f))]
		if not files:
			return 0
		f = files[0]
		t = os.path.getmtime(os.path.join(self.dirname(), f))
		for f in files[1:]:
			t = min(t, os.path.getmtime(os.path.join(self.dirname(), f)))
		return t

	def output_files_exist(self):
		for f in self.output_files:
			path = os.path.join(self.dirname(), f)
			if os.path.exists(path) and os.path.isfile(path):
				continue
			else:
				return False
		return True

	def dir_exists(self):
		return os.path.exists(self.dirname())

	def can_be_skipped(self, timestamp):
		if self.output_files_exist():
			if not timestamp:						#first task
				return True
			elif self.firstmodified() > timestamp:	#later than timestamp
				return True
		return False

	def must_be_rewritten(self, timestamp):
		if self.firstmodified() < timestamp:
			return True
		return False

	def prepare(self, prev):
		halt(not isinstance(prev, self.prev_class), "Wrong previous step: expected %s, got %s" % (self.prev_class, type(prev)))
		self.prev = prev
		if self.prev:
			self.ecp = prev.ecp
			self.basis = prev.basis

	def run(self):
		drn = self.dirname()
		if not os.path.exists(drn):
			os.makedirs(drn)

	def finalize(self):
		pass

class DownloadStep(Step):

	GRECP_FILENAME = 'grecp.inp'
	ECP_FILENAME = 'ecp.inp'
	BASIS_FILENAME = 'basis.inp'

	def __init__(self, grecpfilepath = None):
		Step.__init__(self)
		self.title = 'Download'
		self.description = 'Downloading GRECP and basis'
		self.output_files = [DownloadStep.ECP_FILENAME, DownloadStep.BASIS_FILENAME]
		self.prev_class = type(None)
		self.grecpfilepath = grecpfilepath

	def run(self):
		Step.run(self)
		if self.grecpfilepath:
			shutil.copy(self.grecpfilepath, os.path.join(self.dirname(), DownloadStep.GRECP_FILENAME))
		else:
			logging.info('Enter element symbol or atomic number')
			while True:
				s = input().lower().strip()
				try:
					n = parse_element_string(s)
					break
				except Exception:
					pass
			save_mos_ecp(n, os.path.join(self.dirname(), DownloadStep.GRECP_FILENAME))
		ecp = read_ecp_mos(os.path.join(self.dirname(), DownloadStep.GRECP_FILENAME))
		save_starting_basis(ecp.n, os.path.join(self.dirname(), DownloadStep.BASIS_FILENAME))
		write_ecp_nw(ecp, os.path.join(self.dirname(), DownloadStep.ECP_FILENAME))

	def finalize(self):
		Step.finalize(self)
		logging.info("Semi-local ECP is written to '%s'" % DownloadStep.ECP_FILENAME)
		logging.info("Basis set to take exponents from is written to '%s'" % DownloadStep.BASIS_FILENAME)
		logging.info("Feel free to modify these files.")
		logging.info("Press ENTER to start next task.")
		input()
		self.ecp = read_ecp_nw(os.path.join(self.dirname(), DownloadStep.ECP_FILENAME))
		self.basis = read_basis_nw(os.path.join(self.dirname(), DownloadStep.BASIS_FILENAME))


class OptimizeStep(Step):
	BASIS_OPT_FILENAME = 'exponents.optimized'
	BASIS_START_FILENAME = 'exponents.start'
	BASIS_RESTART_FILENAME = 'exponents.restart'
	CONSTRAINTS_FILENAME = 'opt_constraints'
	STATE_FILENAME = 'opt_state'

	OVERWRITE = "Overwrite and start from previous step result"
	START = "Use 'exponents.start' file"
	RESTART = "Use 'exponents.restart' file"

	def __init__(self):
		Step.__init__(self)
		self.title = 'Optimize'
		self.description = 'Optimizing core exponents'
		self.input_files = [OptimizeStep.BASIS_START_FILENAME]
		self.output_files = [OptimizeStep.BASIS_OPT_FILENAME]
		self.prev_class = DownloadStep

	def run(self):
		Step.run(self)
		variants = [OptimizeStep.OVERWRITE]
		if os.path.exists(os.path.join(self.dirname(), OptimizeStep.BASIS_START_FILENAME)):
			variants.append(OptimizeStep.START)
		if os.path.exists(os.path.join(self.dirname(), OptimizeStep.BASIS_RESTART_FILENAME)):
			variants.append(OptimizeStep.RESTART)
		variant = OptimizeStep.OVERWRITE
		if len(variants) > 1:
			_, variant = multiselect('Previous exponents input found:', variants)
		es = None
		start_filename = OptimizeStep.BASIS_START_FILENAME
		if variant == OptimizeStep.OVERWRITE:
			desired_maxl = estimate_maxl(self.ecp.n, self.ecp.ncore)
			basis_maxl = self.basis.maxl()
			if basis_maxl > desired_maxl:
				logging.warning('Polarization functions are present.')
				cut = yesno('Cut polarization functions with l > %d?' % desired_maxl)
				if cut:
					self.basis.cut_polarization(desired_maxl)
			es = ExponentSet.from_basis(self.basis)
			select_vc_divide(es, self.ecp.n, self.ecp.ncore, desired_maxl)
			es.write_to_file(os.path.join(self.dirname(), OptimizeStep.BASIS_START_FILENAME))
			logging.info("Exponents are written to '%s'" % OptimizeStep.BASIS_START_FILENAME)
		elif variant == OptimizeStep.START:
			pass
		elif variant == OptimizeStep.RESTART:
			start_filename = OptimizeStep.BASIS_RESTART_FILENAME
		if os.path.exists(os.path.join(self.dirname(), OptimizeStep.CONSTRAINTS_FILENAME)):
			pass
		else:
			cs = ConstraintSettings(self.basis.maxl(), {})
			cs.write_to_file(os.path.join(self.dirname(), OptimizeStep.CONSTRAINTS_FILENAME))
			logging.info("Constraints are written to '%s'" % OptimizeStep.CONSTRAINTS_FILENAME)
		if os.path.exists(os.path.join(self.dirname(), OptimizeStep.STATE_FILENAME)):
			pass
		else:
			state = select_atom_state(self.ecp.n, self.ecp.nval, os.path.join(self.dirname(), OptimizeStep.STATE_FILENAME))
		logging.info("Exponents will be taken from '%s'" % start_filename)
		logging.info("Asterisks(*) mark exponents to optimize.")
		logging.info("Constraints will be taken from '%s'" % OptimizeStep.CONSTRAINTS_FILENAME)
		logging.info("Atom state will be taken from '%s'" % OptimizeStep.STATE_FILENAME)
		logging.info("Feel free to modify these files before start.")
		logging.info("Press ENTER to start optimization.")
		input()
		es = ExponentSet.from_file(os.path.join(self.dirname(), start_filename))
		cs = ConstraintSettings.from_file(os.path.join(self.dirname(), OptimizeStep.CONSTRAINTS_FILENAME))
		state = AtomState.from_file(os.path.join(self.dirname(), OptimizeStep.STATE_FILENAME))
		restart_path = os.path.join(self.dirname(), OptimizeStep.BASIS_RESTART_FILENAME)
		lset = None
		while True:
			es = optimize_exponents(es, self.ecp, state, cs, 500, self.dirname(), OptimizeStep.BASIS_RESTART_FILENAME, lset)
			es, discarded_ls = check_for_discarding(es, self.ecp, state, self.dirname())
			if discarded_ls:
				es.write_to_file(os.path.join(self.dirname(), OptimizeStep.BASIS_RESTART_FILENAME))
				lset = discarded_ls
				print_subtitle('Re-optmizing basis after discarding')
			else:
				break
		es.write_to_file(os.path.join(self.dirname(), OptimizeStep.BASIS_OPT_FILENAME))

	def finalize(self):
		Step.finalize(self)
		logging.info("Optimized exponents are written to '%s'" % OptimizeStep.BASIS_OPT_FILENAME)
		logging.info("Editing this file is possible but not recommended.")
		logging.info("Press ENTER to start next task.")
		input()
		self.exponents = ExponentSet.from_file(os.path.join(self.dirname(), OptimizeStep.BASIS_OPT_FILENAME))
		self.state = AtomState.from_file(os.path.join(self.dirname(), OptimizeStep.STATE_FILENAME))
		self.constraints = ConstraintSettings.from_file(os.path.join(self.dirname(), OptimizeStep.CONSTRAINTS_FILENAME))

class ContractStep(Step):

	BASIS_START_FILENAME = 'basis.start'
	BASIS_RESTART_FILENAME = 'basis.restart'
	BASIS_CONTRACTED_FILENAME = 'basis.contracted'

	OVERWRITE = "Overwrite and start from previous step result"
	START = "Use 'basis.start' file"
	RESTART = "Use 'basis.restart' file"


	def __init__(self):
		Step.__init__(self)
		self.title = 'Contract'
		self.description = 'Contracting core exponents'
		self.input_files = [ContractStep.BASIS_START_FILENAME]
		self.output_files = [ContractStep.BASIS_CONTRACTED_FILENAME]
		self.prev_class = OptimizeStep

	def prepare(self, prev):
		Step.prepare(self, prev)
		self.exponents = prev.exponents
		self.state = prev.state
		self.constraints = prev.constraints

	def run(self):
		Step.run(self)
		variants = ['Segmented contraction', 'General contraction', 'No contraction']
		index, variant = multiselect("Select contraction method", variants)
		if index == 0:
			variants = [ContractStep.OVERWRITE]
			if os.path.exists(os.path.join(self.dirname(), ContractStep.BASIS_START_FILENAME)):
				variants.append(ContractStep.START)
			if os.path.exists(os.path.join(self.dirname(), ContractStep.BASIS_RESTART_FILENAME)):
				variants.append(ContractStep.RESTART)
			variant = ContractStep.OVERWRITE
			if len(variants) > 1:
				_, variant = multiselect('Previous basis input found:', variants)
			es = None
			start_filename = ContractStep.BASIS_START_FILENAME
			if variant == ContractStep.OVERWRITE:
				basis = estimate_segmented_contraction(self.exponents, self.ecp, self.state, self.dirname())
				write_basis_nw(basis, os.path.join(self.dirname(), ContractStep.BASIS_START_FILENAME))
				logging.info("Basis is written to '%s'" % ContractStep.BASIS_START_FILENAME)
			elif variant == ContractStep.START:
				pass
			elif variant == ContractStep.RESTART:
				start_filename = ContractStep.BASIS_RESTART_FILENAME
			logging.info("Basis will be taken from '%s'" % start_filename)
			logging.info("Editing this file is possible but not recommended.")
			logging.info("Press ENTER to start optimization.")
			input()
			self.basis = read_basis_nw(os.path.join(self.dirname(), start_filename))
			self.basis = optimize_contraction(self.basis, self.ecp, self.state, 500, self.dirname(), ContractStep.BASIS_RESTART_FILENAME)
			write_basis_nw(self.basis, os.path.join(self.dirname(), ContractStep.BASIS_CONTRACTED_FILENAME))
		elif index == 1:
			self.basis = estimate_general_contraction(self.exponents, self.ecp, self.state, self.dirname())
			write_basis_nw(self.basis, os.path.join(self.dirname(), ContractStep.BASIS_CONTRACTED_FILENAME))
		else:
			self.basis =  self.exponents.build_basis()
			write_basis_nw(self.basis, os.path.join(self.dirname(), ContractStep.BASIS_CONTRACTED_FILENAME))


	def finalize(self):
		Step.finalize(self)
		self.basis = read_basis_nw(os.path.join(self.dirname(), ContractStep.BASIS_CONTRACTED_FILENAME))

class FinalStep(Step):

	NWCHEM_FILENAME = 'basis.nw'
	PYSCF_FILENAME = 'basis.pyscf'
	CRYSTAL_FILENAME = 'basis.crystal'


	def __init__(self):
		Step.__init__(self)
		self.title = 'Output'
		self.description = 'Writing basis and ECP in different formats'
		self.output_files = [FinalStep.NWCHEM_FILENAME, FinalStep.PYSCF_FILENAME, FinalStep.CRYSTAL_FILENAME]
		self.prev_class = ContractStep

	def prepare(self, prev):
		Step.prepare(self, prev)

	def run(self):
		Step.run(self)
		steal_polarization(self.ecp, self.basis)
		nw_ecp = ecp_to_nw(self.ecp)
		py_ecp = ecp_to_pyscf(self.ecp)
		nw_basis = basis_to_nw(self.basis)
		with open(os.path.join(self.dirname(), FinalStep.NWCHEM_FILENAME), 'w') as f:
			f.write(nw_ecp)
			f.write('\n')
			f.write(nw_basis)
		with open(os.path.join(self.dirname(), FinalStep.PYSCF_FILENAME), 'w') as f:
			f.write(py_ecp)
			f.write('\n')
			f.write(nw_basis)
		estimate_and_write_crystal(self.ecp, self.basis, os.path.join(self.dirname(), FinalStep.CRYSTAL_FILENAME))


        #=======================================#
#=======#                  MAIN                 #=======#
        #=======================================#


filename = None
if len(sys.argv) > 1:
	if sys.argv[1][0] != '-':
		filename = sys.argv[1]

sm = StepManager([
		DownloadStep(filename),
		OptimizeStep(),
		ContractStep(),
		FinalStep()
	])

sm.run()


















