import re, sys, math

from ..atom.shells import Shells
from .gauss_function import GaussFunctionContracted, GaussFunctionNormed, GaussFunction, GaussBasis, GaussPseudoPotential


class GaussFormat:

	@staticmethod
	def basis_to_gaussian94(components):
		res = ''
		for l in sorted(components.keys()):
			for cg in components[l]:
				res += "%s  %2d 1.00   \n" % (Shells.SHELLS[l], len(cg.fs))
				for n, gf in cg.fs:
					res += "    %25.18f    %25.18f\n" % (gf.a, n)
		return res

	@staticmethod
	def basis_from_gaussian94(res):
		basis = GaussBasis()
		gs = []
		cg = None
		l = -1
		for line in res[1:]:
			#TODO: REGEXP
			if line[0] != ' ' and line[0] != '	':
					if cg:
						basis.add_function(cg)
					cg = GaussFunctionContracted()
					l = Shells.SHELLS.index(line[0])
			else:
				ls = re.split("\s*", line.strip())
				a = float(ls[0])
				n = float(ls[1])
				cg.fs.append((n, GaussFunctionNormed(a, l)))
		if cg:
			basis.add_function(cg)
		return gs


















