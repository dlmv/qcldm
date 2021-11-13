import re, os, logging

from ..atom.shells import Shells
from ..gauss_functions.gauss_function import GaussFunctionNormed, GaussFunctionContracted

from turbo_format import TurboTemplate

class TurboBasis:

	class TurboEcp:
		def __init__(self, name, ncore):
			self.name = name
			self.ncore = ncore

	def __init__(self, name, functions):
		self.name = name
		self.functions = functions

	def orbnum(self):
		res = 0
		for l in range(len(Shells.SHELLS)):
			for f in self.functions:
				if f.fs[0][1].l == l:
					res += 2 * l + 1
		return res

	def orbarray(self):
		res = []
		for l in range(len(Shells.SHELLS)):
			nl = 0
			for f in self.functions:
				if f.fs[0][1].l == l:
					nl += 2 * l + 1
			res.append(nl)
		return res

	def internal_overlap(self):
		olp = []
		for f1 in self.functions:
			line = []
			for f2 in self.functions:
				line.append(f1.overlap(f2))
			olp.append(line)
		return olp

	@staticmethod
	def read_basis(name):
		tt = TurboTemplate.from_file(name)
		bases = {}
		n = 0
		lines = tt.param("basis").multiparam
		pattern = re.compile("\s*([0-9]+)\s*([spdfgh])")
		while n < len(lines):
			l = lines[n]
			if l.strip() == '*':
				n += 1
				if n >= len(lines):
					break
				name = lines[n].strip()
				n += 1
				basis_def = lines[n].split()
				ls = basis_def[-1][1:-1].split('/')
				n += 2
				functions = []
				while n < len(lines) and lines[n].strip() != '*':
					gc = GaussFunctionContracted()
					m = pattern.match(lines[n])
					if m:
						ng = int(m.group(1))
						l = Shells.SHELLS.lower().index(m.group(2))
						for n in range(n + 1, n + ng + 1):
							a, c = [float(x) for x in lines[n].split()]
							g = GaussFunctionNormed(a, l)
							gc.fs.append([c, g])
						gc.normalize()
						functions.append(gc)
					n += 1
				tb = TurboBasis(name, functions)
				bases[name] = tb
		
		ecps = {}
		lines = tt.param("ecp").multiparam
		n = 0
		while n < len(lines):
			l = lines[n]
			if l.strip() == '*':
				n += 1
				if n >= len(lines):
					break
				name = lines[n].strip()
				n += 2
				ecp_def = lines[n].strip()
				ecp_def_map = {}
				for m in re.finditer('(\\w+)=\s+(\d+)', ecp_def):
					ecp_def_map[m.group(1)] = int(m.group(2))
				te = TurboBasis.TurboEcp(name, ecp_def_map['ncore'])
				ecps[name] = te

				
				while n < len(lines) and lines[n].strip() != '*':
					n += 1
				if 'lsomax' in ecp_def_map.keys():
					n += 1
					while n < len(lines) and lines[n].strip() != '*':
						n += 1
		return bases, ecps


















