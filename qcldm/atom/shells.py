from ..util.elements import ELEMENTS

class Shells:
	SHELLS = "SPDFGH"
	SHELL_POP = [(4 * l + 2) for l in range(len(SHELLS))]

	@staticmethod
	def aufbau_sequence():
		nl = 1
		while True:
			for l in reversed(range(nl)):
				n = nl - l
				if n > l:
					yield n, l
			nl += 1

	@staticmethod
	def direct_sequence():
		n = 1
		while True:
			for l in range(n):
				yield n, l
			n += 1

	@staticmethod
	def estimate_valence(z):
		eld = ELEMENTS[z].eleconfig_dict
		nmax = 0
		for n, ls in eld.keys():
			if n > nmax:
				nmax = n
		val = 0
		for n, ls in eld.keys():
			l = Shells.SHELLS.lower().index(ls)
			if n + max(0, l - 1) == nmax:
				val += eld[n, ls]
		return val

	@staticmethod
	def estimate_pf_pop(z, zval):
		zcore = z - zval
		zcur = 0
		cores = {}
		for n, l in Shells.direct_sequence():
			if l not in cores.keys():
				cores[l] = 0
			cores[l] += 1
			zcur += Shells.SHELL_POP[l]
			if zcur >= zcore:
				break
		assert zcur == zcore, "Cannot get core for %d electrons" % zcore
		eld = ELEMENTS[z].eleconfig_dict
		pops = {}
		for n, ls in eld.keys():
			l = Shells.SHELLS.lower().index(ls)
			nn = n
			if l in cores.keys():
				nn -= cores[l]
			if nn > l:
				pops[(nn, l)] = eld[(n, ls)]
		return pops
