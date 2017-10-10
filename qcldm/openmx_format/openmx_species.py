import re, os, logging
from ..util.units import Units
from ..atom.shells import Shells
from pao_format import PAO
from vps_format import VPS

class openmx_species:

	def __init__(self, ls, data_path):
		self.numorbs = {}
		lsp = filter(None, re.split('-', ls[1]))
		lso =  re.findall("[a-z]\d*", lsp[1])
		for o in lso:
			oname = o[0]
			onumber = int(o[1:])
			self.numorbs[oname] = onumber
		self.name = ls[0]
		self.pao = lsp[0]
		self.vps = ls[2]

		lsr = re.findall("\d*\.\d*", lsp[0])
		self.paor = float(lsr[0])

		pao_path = os.path.join(data_path, "PAO", self.pao+".pao")
		vps_path = os.path.join(data_path, "VPS", self.vps+".vps")
		self.pp = None
		self.basis = None
		try:
			self.pp = VPS.from_file(vps_path)
			self.basis = PAO.from_file(pao_path)
		except IOError:
			logging.warn(u'PAO or VPS file not found for {}'.format(self.name))

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

	def fulllist(self, prefix=''):
		orbstr = ''
		for s in Shells.SHELLS.lower():
			if s in self.numorbs.keys():
				orbstr += s + str(self.numorbs[s])
		return [prefix + self.name,  prefix + self.pao + '-' + orbstr, prefix + self.vps]

	def prefixed(self, prefix, path):
		return openmx_species(self.fulllist(prefix), path)

	def real_r(self):
		return self.paor * Units.BOHR / Units.UNIT

	def estimate_valence(self):
		return Shells.estimate_valence(self.pp.z)

