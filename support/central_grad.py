#!/usr/bin/python 
import os, sys

def read_embed_positions():
	pos = []
	nepos = []
	with open("coord") as c:
		n = 0
		for l in c:
			ls = l.split()
			if '$' not in l and ls:
				if ls[-1] == 'zz':
					pos += [n]
				elif ls[-1] == 'ne':#FIXME
					nepos += [n]
				n += 1
	return n, pos, nepos

def read_grad_from_control(num, pos, nepos):
	with open("control") as c:
		lines = c.read().splitlines()
		n = 0
		while n < len(lines):
			if '$grad' in lines[n]:
				break
			n += 1
		assert n < len(lines)
		grad_s = 0.
		grad_a = 0.
		while n < len(lines):
			while n < len(lines):
				if lines[n][0] == '#':
					n += 1
					continue
				if 'cycle' in lines[n]:
					break
				n += 1
			if n == len(lines):
				return grad_s, grad_a
			grad_s = 0.
			grad_a = 0.
			rawgrads = []
			n0 = n
			for i in range(num):
				n = n0 + num + 1 +  i
				igrads = [float(l.replace("D", "E")) for l in lines[n].split()]
				if i not in pos + nepos:
					rawgrads.extend(igrads)
			grad_s = (reduce(lambda s, i: s + i**2, rawgrads, 0.))**0.5 
			grad_a = grad_s / (len(rawgrads) / 3)**0.5


n, pos, nepos = read_embed_positions()
grad_s, grad_a = read_grad_from_control(n, pos, nepos)
print grad_s
print grad_a
























