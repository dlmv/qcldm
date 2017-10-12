import re, os, logging
from ..atom.shells import Shells

from turbo_format import TurboTemplate

class MosReader:
	def __init__(self):
		self.base_format = None
		self.matrix = []

	def load(self, of):
		self.base_format = of
		assert self.base_format.data[0].key.startswith('uhfmo')
		fmt = self.base_format.data[0].lineparam.split()[-1]
		m = re.match('format\(([0-9]+)d([0-9]+)\.([0-9]+)\)', fmt)
		assert m
		nd = int(m.group(1))
		dlen = int(m.group(2))
		lines = self.base_format.data[0].multiparam[4:]
		n = 0
		self.matrix = []
		while n < len(lines):
			mpart = []
			l = lines[n]
			nsaos = int(l.split()[-1].split('=')[-1])
			ns = nsaos / nd
			for n in range(n + 1, n + ns + 1):
				for i in range(len(lines[n]) / dlen):
					d = lines[n][i * dlen:(i + 1) * dlen].replace('D', 'E')
					f = float(d)
					mpart.append(f)
			self.matrix.append(mpart)
			n += 1
		print self.matrix


	@staticmethod
	def from_string(datastring):
		res = MosReader()
		res.load(TurboTemplate.from_string(datastring))
		return res

	@staticmethod
	def from_file(name):
		res = MosReader()
		res.load(TurboTemplate.from_file(name))
		return res
