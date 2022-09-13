from ..util.elements import ELEMENTS

ANIONS = [7, 8, 9, 16, 17, 35, 53]#FIXME: this is ugly :/

class Shells:
	SHELLS = "SPDFGH"
	SHELL_POP = [(4 * l + 2) for l in range(len(SHELLS))]

	@staticmethod
	def aufbau_sequence():
		nl = 1
		while True:
			for l in reversed(list(range(nl))):
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
		for n, ls in list(eld.keys()):
			if n > nmax:
				nmax = n
		val = 0
		for n, ls in list(eld.keys()):
			l = Shells.SHELLS.lower().index(ls)
			if n + max(0, l - 1) == nmax:
				val += eld[n, ls]
		return val

	@staticmethod
	def estimate_valence_byname(a):
		try:
			z = ELEMENTS[a].number
			return Shells.estimate_valence(z)
		except Exception:
			return 0

	@staticmethod
	def estimate_charge_byname(a):
		z = ELEMENTS[a].number
		return Shells.estimate_charge(z)

	@staticmethod
	def estimate_charge(z):
		if z in ANIONS:
			return Shells.estimate_valence(z) - 8
		else:
			return Shells.estimate_valence(z)

	@staticmethod
	def estimate_pf_pop(z, zval):
		zcore = z - zval
		zcur = 0
		cores = {}
		for n, l in Shells.direct_sequence():
			if l not in list(cores.keys()):
				cores[l] = 0
			cores[l] += 1
			zcur += Shells.SHELL_POP[l]
			if zcur >= zcore:
				break
		assert zcur == zcore, "Cannot get core for %d electrons" % zcore
		eld = ELEMENTS[z].eleconfig_dict
		pops = {}
		for n, ls in list(eld.keys()):
			l = Shells.SHELLS.lower().index(ls)
			nn = n
			if l in list(cores.keys()):
				nn -= cores[l]
			if nn > l:
				pops[(nn, l)] = eld[(n, ls)]
		return pops
