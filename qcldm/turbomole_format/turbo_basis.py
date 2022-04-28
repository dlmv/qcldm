import re, os, logging

from ..atom.shells import Shells
from ..gauss_functions.gauss_function import GaussFunction, GaussFunctionNormed, GaussFunctionContracted

from turbo_format import TurboTemplate

class TurboBasis:

	class EcpPart:
		def __init__(self, functions):
			self.functions = functions

	class TurboEcp:
		def __init__(self, name, ncore, local, semilocal, spinorbit):
			self.name = name
			self.ncore = ncore
			self.local = local
			self.semilocal = semilocal
			self.spinorbit = spinorbit

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
				while '#' == lines[n][0]:
					n += 1
				n += 1
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
				n += 1
				while '#' in lines[n]:
					n += 1
				n += 1
				ecp_def = lines[n].strip()
				ecp_def_map = {}
				for m in re.finditer('(\\w+)=\s+(\d+)', ecp_def):
					ecp_def_map[m.group(1)] = int(m.group(2))

				n += 1
				while '#' == lines[n][0]:
					n += 1
				lmax = ecp_def_map['lmax']
				assert lines[n].strip().split()[0] == Shells.SHELLS[lmax].lower(), 'wrong lmax for %s' % name
				n += 1
				tmpfuncs = []
				while len(lines[n].split()) == 3 and lines[n].split()[1] in ['0', '1', '2']:
					ls = lines[n].split()
					c, l, a = float(ls[0]), int(ls[1]), float(ls[2])
					tmpfuncs.append([c, GaussFunction(a, l)])
					n += 1
				local = TurboBasis.EcpPart(tmpfuncs)
				semilocal = []
				
				for l in range(ecp_def_map['lmax']):
					assert lines[n].strip().split()[0] == '%s-%s' % (Shells.SHELLS[l].lower(), Shells.SHELLS[lmax].lower()), 'wrong l for %s' % name
					n += 1
					tmpfuncs = []
					while len(lines[n].split()) == 3 and lines[n].split()[1] in ['0', '1', '2']:
						ls = lines[n].split()
						c, l, a = float(ls[0]), int(ls[1]), float(ls[2])
						tmpfuncs.append([c, GaussFunction(a, l)])
						n += 1
					semilocal.append(TurboBasis.EcpPart(tmpfuncs))
				assert lines[n].strip() == '*'
				spinorbit = []
				if 'lsomax' in ecp_def_map.keys():
					n += 1
					for l in range(1, ecp_def_map['lsomax'] + 1):
						assert lines[n].strip().split()[0] == '%s-SpinOrbit' % Shells.SHELLS[l].lower(), 'wrong so for %s' % name
						n += 1
						tmpfuncs = []
						while len(lines[n].split()) == 3 and lines[n].split()[1] in ['0', '1', '2']:
							ls = lines[n].split()
							c, l, a = float(ls[0]), int(ls[1]), float(ls[2])
							tmpfuncs.append([c, GaussFunction(a, l)])
							n += 1
						spinorbit.append(TurboBasis.EcpPart(tmpfuncs))	
				assert lines[n].strip() == '*'
					
				te = TurboBasis.TurboEcp(name, ecp_def_map['ncore'], local, semilocal, spinorbit)
				ecps[name] = te
		return bases, ecps


















