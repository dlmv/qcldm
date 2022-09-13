import math, logging, sys
from ..util.elements import ELEMENTS
from .numeric_function import NumericFunction, NumericOperations
from shells import Shells


class NumericBasis:

	def __init__(self):
		self.z = 0
		self.etotal = 0
		self.eval = 0
		self.cutoff = 0
		self.functions = {}
		self.vden = None
		self.grid = []

	def load_ze_from_pp(self, pp):
		self.z = pp.z
		self.etotal = pp.etotal
		self.eval = pp.eval

	def integrated_vden(self):
		return NumericOperations.integrate(self.vden.data) * math.pi * 4

	def estimate_valence(self):
		return Shells.estimate_valence(self.z)
	
	def estimate_full_vden(self): 
		vden = 0 * self.functions[0][1]
		pops = Shells.estimate_pf_pop(self.z, self.eval)
		for n, l in list(pops.keys()):
			f = self.functions[l][n]
			occ = pops[(n,l)]
			vden += occ * (f*f) / math.pi / 4
		return vden

	def estimate_upper_vden(self):
		vden = 0 * self.functions[0][1]
		pops = Shells.estimate_pf_pop(self.z, Shells.estimate_valence(self.z))
		for n, l in list(pops.keys()):
			f = self.functions[l][n]
			occ = pops[(n,l)]
			vden += occ * (f*f) / math.pi / 4
		return vden

	def set_functions(self, fs, grid=None):
		if not grid:
			grid = NumericOperations.loggrid(fs[0].grid()[0], fs[0].grid()[-1])
		self.grid = grid
		self.cutoff = grid[-1]
		self.functions = NumericOperations.resort_ln([abs(f.rescaled(self.grid)) for f in fs])
		self.shift_functions()

	def shift_functions(self):
		new_functions = {}
		for l in list(self.functions.keys()):
			new_functions[l] = {}
			ns = sorted(self.functions[l].keys())
			nn = l + 1
			for n in ns:
				new_functions[l][nn] = self.functions[l][n]
				nn += 1
		self.functions = new_functions

	def ortonorm(self):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Ortonorming basis functions')
		logging.info('*********************************************')
		logging.info('')
		for l in sorted(self.functions.keys()):
			logging.debug('L=%d' % l)
			fs = []
			for n in sorted(self.functions[l].keys()):
				fs.append(self.functions[l][n])
			self.functions[l] = {}
			fs = NumericOperations.gramm_schmidt(fs)
			for n, f in enumerate(fs):
				if f.norm() > 0.05:
					f.normalize()
					self.functions[l][f.n] = f
				else:
					logging.warn(" %d'th basis function for L=%d is linearly dependent" % (n, l))
		self.shift_functions()

	def align_functions(self, tol=0.7):
		n0 = len(self.functions[0])
		for l in sorted(self.functions.keys()):
			if 1. * len(self.functions[l]) / n0 < tol:
				break
		lmax = l
		n = sys.maxsize
		for l in sorted(self.functions.keys()):
			if l <= lmax:
				n = min(n, len(self.functions[l]))
		nmax = n
		new_functions = {}
		for l in sorted(self.functions.keys())[:lmax + 1]:
			new_functions[l] = {}
			for n in sorted(self.functions[l].keys())[:nmax]:
				new_functions[l][n] = self.functions[l][n]
		self.functions = new_functions










