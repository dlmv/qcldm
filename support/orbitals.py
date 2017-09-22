#!/usr/bin/python
import os, re, sys

start = 'Real (Re) and imaginary (Im) parts of LCAO coefficients'
homo = '  HOMO = '

class AO:
	def __init__(self, atom, shell):
		self.atom = atom #(1, H)
		self.shell = shell #(0, s)

class LCAO:
	def __init__(self, n, e):
		self.n = n
		self.e = e
		self.aos = []

	def aotype(self):
		ls = sorted(self.aos, key=lambda tup: tup[0], reverse=True)
		res = []
		for l in ls:
			if l[0] >= ls[0][0] * 0.5:
				res.append(l[1])
		return res

	def aotypestring(self):
		res = ''
		for l in self.aotype():
			res += str(l.atom[0]) + l.atom[1] + ':' + str(l.shell[0]) + l.shell[1] + ', '
#			res += str(l.shell[0]) + l.shell[1] + ', '
		return res[:-2]

	def aosimplestring(self):
		l = self.aotype()[0]
		res = str(l.shell[0]) + l.shell[1][0]
		return res

def find_start(lines, n):
	while n < len(lines):
		if start in lines[n]:
			return n
		n += 1
	return 0

def find_homo(lines, n):
	while n < len(lines):
		if homo in lines[n]:
			return n
		n += 1
	return 0

def read_line(line, atom, hnum):
	ls = filter(None, re.split('\s*', line[:-1]))
	if len(ls) - hnum * 4 == 4:
		atom = (int(ls[0]), ls[1])
		ls = ls[2:]
	elif len(ls) - hnum * 4 != 2:
		raise Exception('Strange line:\n %s' % line)
	shell = (int(ls[0]), ls[1])
	ls = ls[2:]
	o = AO(atom, shell)
	caos = []
	for n in range(hnum):
		k = 0
		for c in ls[n*4:(n+1)*4]:
			k += float(c)**2
		caos.append((k, o))
	return caos, atom

def read_portion(lines, n):
	found = False
	lcaos = []
	while n < len(lines) and not found:
		if len(lines[n]) < 3:
			n += 1
		elif '*****' not in lines[n]:
			found = True
		else:
			return None, len(lines)
	if found:
		nums = filter(None, re.split('\s*', lines[n][:-1]))
		n += 1
		es = filter(None, re.split('\s*', lines[n][:-1]))
		n += 2
		cs = filter(None, re.split('\s*', lines[n][:-1]))
		n += 2
		for k in range(len(nums)):
			lcaos.append(LCAO(int(nums[k]), float(es[k])))
		stop = False
		atom = None
		while n < len(lines) and not stop:
			caos, atom = read_line(lines[n], atom, len(nums))
			for l, cao in zip(lcaos, caos):
				l.aos.append(cao)
			n += 1
			if len(lines[n]) < 3:
				stop = True
	return lcaos, n

def print_string(l, n, nocc, verbose):
	occ = '*' if n < nocc else ' '
	return '%s%3d %9.5f %s' % (occ, l.n, l.e, l.aotypestring() if verbose else l.aosimplestring())

def parse_orbitals(name, verbose=False):
	f = open(name)
	lines = f.readlines()
	f.close()
	n = find_homo(lines, 0)
	nocc = int(lines[n][len(homo):])
	n = find_start(lines, n)
	lcaos = []
	while n < len(lines):
		tlcaos, n = read_portion(lines, n + 1)
		if tlcaos:
			lcaos.extend(tlcaos)
	for n in range(len(lcaos)):
		l = lcaos[n]
		s = print_string(l, n, nocc, verbose)
		if n == 0 or l.e != lcaos[n-1].e:
			print s

verbose = '-v' in sys.argv[2:]
parse_orbitals(sys.argv[1], verbose)

