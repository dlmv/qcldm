#!/usr/bin/python
import sys, os, re
sys.dont_write_bytecode = True

float_regex = "[-+]?[0-9]*\.?[0-9]+"

shells = "SPDFGH"

def parse_ljmj(s):
	m = re.match("\|([SPDFGH]),(%s),(%s)\>" % (float_regex, float_regex), s)
	l = shells.index(m.group(1))
	j = float(m.group(2))
	mj = float(m.group(2))
	return l, j, mj

def read_data(root):
	dm = []
	for dr, subFolders, files in os.walk(root):
		dname = os.path.basename(os.path.normpath(dr))
		if dname.startswith('JReduced'):
			if dm:
				print 'Multiple logs in %s' % root
				print 'Aborting!'
				return
			with open(os.path.join(os.path.normpath(dr), 'JRed_DM_Re_notAcc_')) as f:
				lines = f.read().splitlines()
				N = int(lines[0])
				seq = []
				header = lines[1].split()
				ljs = {}
				ls = set()
				for h in header:
					l, j, mj = parse_ljmj(h)
					seq.append((l, j, mj))
					ljs[(l, j)] = 0
					ls.add(l)
				for n in range(N):
					d = float(lines[n + 2].split()[n + 1])
					l, j, mj = seq[n]
					ljs[(l, j)] += d
				for l in sorted(ls):
					if l != 0:
						print l, ljs[(l, l + 0.5)] / ljs[(l, l - 0.5)]
#				for l, j in sorted(ljs.keys()):
#					print l, j, ljs[(l, j)]


root = sys.argv[1]

read_data(root)



















