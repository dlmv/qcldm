import re
from ..atom.shells import Shells
from ..functions.gauss_function import GaussFunctionContracted, GaussFunctionNormed


class GaussFormat:

	@staticmethod
	def to_gaussian94(gs):
		res = ''
		for cg in gs:
			res += "%s  %2d 1.00   \n" % (Shells.SHELLS[cg.fs[0][1].l], len(cg.fs))
			for n, gf in cg.fs:
				res += "    %25.18f    %25.18f\n" % (gf.a, n)
		return res

	@staticmethod
	def from_gaussian94(res):
		gs = []
		cg = None
		l = -1
		for line in res[1:]:
			#TODO: REGEXP
			if line[0] != ' ' and line[0] != '	':
					if cg:
						gs.append(cg)
					cg = GaussFunctionContracted()
					l = Shells.SHELLS.index(line[0])
			else:
				ls = re.split("\s*", line.strip())
				a = float(ls[0])
				n = float(ls[1])
				cg.fs.append((n, GaussFunctionNormed(a, l)))
		if cg:
			gs.append(cg)
		return gs
