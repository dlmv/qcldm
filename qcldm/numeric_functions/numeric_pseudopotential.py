import logging
from ..util.elements import ELEMENTS
from .numeric_function import NumericFunction, NumericOperations, CustomFunctions
from .numeric_potential import create_pseudo_potential, blochl_form


class SeparablePseudoPotential:

	def __init__(self):
		self.z = 0
		self.etotal = 0
		self.eval = 0
		self.loc = NumericFunction()
		self.blochl_projectors = []
		self.pcc = None
		self.grid = []

	def name(self):
		return ELEMENTS[self.z].symbol

class GRECPPseudoPotential:

	def __init__(self, pps, pfs, loc, z, zval):
		self.z = z
		self.zval = zval
		self.loc = loc
		self.grid = loc.grid()
		self.pps = NumericOperations.resort_lnj([f.rescaled(self.grid) for f in pps])
		self.pfs = NumericOperations.resort_lnj([abs(f.rescaled(self.grid)) for f in pfs])
		self.projs = {}
		self.PPS = {}

	def set_projectors(self, projs):
		self.projs = NumericOperations.resort_lnj([f.rescaled(self.grid) for f in projs])

	def name(self):
		return ELEMENTS[self.z].symbol

	def convert_potentials(self, grecp):
		form = "GRECP" if grecp else "semilocal"
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Converting raw PPs to %s form' % form)
		logging.info('*********************************************')
		self.PPS = {}
		for l in list(self.pps.keys()):
			logging.debug('L=%d' % l)
			self.PPS[l] = {}
			nlist = sorted(self.pps[l].keys())
			nv = nlist[-1]
			logging.debug(' valent N=%d' % nv)
			v0 = self.pps[l][nv][l + 0.5]
			v1 = self.pps[l][nv][l - 0.5] if l != 0 else v0

			c0 = []
			c1 = []
			cs = ""
			for nc in nlist[:-1]:
				c0.append([self.pps[l][nc][l + 0.5], self.pfs[l][nc][l + 0.5]])
				if l != 0:
					c1.append([self.pps[l][nc][l - 0.5], self.pfs[l][nc][l - 0.5]])
				cs += (str(nc) + ", ")
			cs = cs[:-2]
			logging.debug(' GRECP core shells: %s' % (cs if cs else 'none'))
			self.PPS[l][l + 0.5] = create_pseudo_potential(v0, c0, self.loc, grecp)
			if l != 0:
				self.PPS[l][l - 0.5] = create_pseudo_potential(v1, c1, self.loc, grecp)

	def create_separable_form(self, blochl_proj, sine_proj=0):
		blochl_p = self.prepare_blochl_form(blochl_proj, sine_proj)
		spp = SeparablePseudoPotential()
		spp.z =self.z
		spp.etotal = self.z
		spp.eval = self.zval
		spp.grid = self.grid
		spp.loc = self.loc
		spp.blochl_projectors = blochl_p
		spp.pcc = None
		return spp

	def expand_vf(self, p, fs, n):
		nfs = []
		for f in fs:
			tmp = f
			for i in range(n):
				nfs.append(tmp)
				tmp = p(tmp)
		return nfs

	def prepare_blochl_form(self, blochl_proj, sine_proj):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Creating Blochl Projectors')
		logging.info('*********************************************')
		logging.info('')

		projectors = []
		proj_energies = []
		proj_ls = []

		for l in sorted(self.pps.keys()):
			logging.debug('L=%d' % l)
			p = self.PPS[l][l + 0.5], self.PPS[l][(l - 0.5) if l != 0 else (l + 0.5)]
			fs = [[],[]]
			nv = sorted(self.pps[l].keys())[-1]
			for n in sorted(self.pps[l].keys()):
				if l in list(self.pfs.keys()) and n in list(self.pfs[l].keys()):
					logging.debug(' Adding pseudofunction for N=%d' % n)
					f0 = self.pfs[l][n][l + 0.5]
					f1 = self.pfs[l][n][l - 0.5] if l != 0 else f0
					fs[0].append(f0)
					fs[1].append(f1)

			if blochl_proj > 1:
				logging.debug(' Expanding projector basis by factor of %d' % blochl_proj)
				fs[0] = self.expand_vf(p[0], fs[0], blochl_proj)
				fs[1] = self.expand_vf(p[1], fs[1], blochl_proj)


			if l in list(self.projs.keys()):
				logging.debug(' Adding %d explicit projectors' % len(self.projs[l]))
				for n in sorted(self.projs[l].keys()):
					f0 = self.projs[l][n][l + 0.5]
					f1 = self.projs[l][n][l - 0.5] if l != 0 else f0
					fs[0].append(f0)
					fs[1].append(f1)

			pmax = 0
			for f in fs[0]:
				pmax = max(pmax, p[0](f).get_cutoff(0.001))
			for f in fs[1]:
				pmax = max(pmax, p[1](f).get_cutoff(0.001))

			if sine_proj > 0:
				logging.debug(' Adding %d sine waves in R=%f' % (sine_proj, pmax))
				for nn in range(sine_proj):
					f0 = CustomFunctions.restricted_sine(self.grid, pmax, nn)
					f1 = CustomFunctions.restricted_sine(self.grid, pmax, nn)
					fs[0].append(f0)
					fs[1].append(f1)

			tp2 = [[],[]]
			te2 = [[],[]]
			for j in range(2):
				tp, te = blochl_form(l, fs[j], p[j], False)
				tp2[j] = tp
				te2[j] = te
			for p1, p2, e1, e2 in zip(tp2[0], tp2[1], te2[0], te2[1]):
				projectors.append((p1, p2))
				proj_energies.append((e1, e2))
				proj_ls.append(l)
		return list(zip(proj_ls, projectors, proj_energies))


	




