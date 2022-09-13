#!/usr/bin/python

import sys, os, re
sys.dont_write_bytecode = True

shells = 'spdfgh'

def orbcompare(o1, o2):
	if o1.n != o2.n:
		return -1 if o1.n < o2.n else 1
	if o1.l != o2.l:
		return -1 if o1.l < o2.l else 1
	if o1.j2 != o2.j2:
		return -1 if o1.j2 < o2.j2 else 1
	return 0

class Orbital:
	def __init__(self):
		self.n = 0
		self.l = 0
		self.j2 = 0
		self.e = 0.

	def __str__(self):
		return "%d%s%d/2" % (self.n, shells[self.l], self.j2)

	def __repr__(self):
		return self.__str__()

	def __lt__(self, other):
		return orbcompare(self, other) < 0
	def __gt__(self, other):
		return orbcompare(self, other) > 0
	def __eq__(self, other):
		return orbcompare(self, other) == 0
	def __le__(self, other):
		return orbcompare(self, other) <= 0
	def __ge__(self, other):
		return orbcompare(self, other) >= 0
	def __ne__(self, other):
		return orbcompare(self, other) != 0

def selection_rule(o1, o2):
	if o1.n >= o2.n:
		return False
	if abs(o1.l - o2.l) != 1:
		return False
	if abs(o1.j2 - o2.j2) != 2 and abs(o1.j2 - o2.j2) != 0:
		return False
	return True


class Transition:
	def __init__(self):
		self.o1 = None
		self.o2 = None
		self.e = 0.

	def __eq__(self, other):
		return self.o1 == other.o1 and self.o2 == other.o2
	def __ne__(self, other):
		return not self.__eq__(other)

	def __str__(self):
		return "%s->%s" % (self.o2, self.o1)

	def __repr__(self):
		return self.__str__()
	

file_regex = 'JRed_Prop_ljm electron part of ([0-9])([%s])([0-9])_([0-9]) orbital energy shift' % shells
energy_regex = '\s*value Re \(ljm electron part of .*? orbital energy shift\) = ([.0-9]*)'



def read_data(root):
	orbitals = []
	found = False
	for dr, subFolders, files in os.walk(root):
		dname = os.path.basename(os.path.normpath(dr))
		if dname.startswith('JReduced'):
			found = True
			if orbitals:
				print('Multiple logs in %s' % root)
				print('Aborting!')
				return
			for f in files:
				r = re.match(file_regex, f)
				if r:
					o = Orbital()
					o.n = int(r.groups()[0])
					o.l = shells.index(r.groups()[1])
					o.j2 = int(r.groups()[2])
					orbitals.append(o)
					with open(os.path.join(dr, f)) as ff:
						for line in ff.read().splitlines():
							rr = re.match(energy_regex, line)
							if rr:
								o.e = float(rr.groups()[0])
	if not found:
		print('No logs in %s' % root)
		print('Aborting!')
		return
	orbitals = sorted(orbitals)
	transitions = []
	for o1 in orbitals:
		for o2 in orbitals:
			if selection_rule(o1, o2):
				t = Transition()
				t.o1 = o1
				t.o2 = o2
				t.e = o2.e - o1.e
				transitions.append(t)
	return transitions

root1 = sys.argv[1]
root2 = sys.argv[2]
trs1 = read_data(root1)
trs2 = read_data(root2)

#for t in trs1:
#	print t, t.e

if trs1 and trs2 and trs1 == trs2:
	for t1, t2 in zip(trs1, trs2):
#		pass
		print(t1, (t2.e - t1.e) * 27211)
























