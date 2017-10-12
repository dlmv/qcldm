import re, os, logging
from ..atom.shells import Shells

from turbo_format import TurboTemplate

class TurboBasis:

	def __init__(self, name, numorbs):
		self.numorbs = numorbs

	def orbnum(self):
		res = 0
		for s in Shells.SHELLS.lower():
			if s in self.numorbs.keys():
				res += self.numorbs[s] * (2 * Shells.SHELLS.lower().index(s) + 1)
		return res

	def orbarray(self):
		res = []
		for s in Shells.SHELLS.lower():
			if s in self.numorbs.keys():
				res.append(self.numorbs[s] * (2 * Shells.SHELLS.lower().index(s) + 1))
			else:
				res.append(0)
		return res

	@staticmethod
	def read_basis(name):
		tt = TurboTemplate.from_file(name)
		res = {}
		n = 0
		lines = tt.param("basis").multiparam
		while n < len(lines):
			l = lines[n]
			if l.strip() == '*':
				n += 1
				name = lines[n].strip()
				n += 1
				basis_def = lines[n].split()
				ls = basis_def[-1][1:-1].split('/')
				numorbs = {}
				for n, x in enumerate(ls):
					numorbs[Shells.SHELLS[n].lower()] = len(x)
				tb = TurboBasis(name, numorbs)
				res[name] = tb
				n += 2
				while n < len(lines) and lines[n].strip() != '*':
					n += 1
			n += 1
		return res
