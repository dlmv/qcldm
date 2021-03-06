import re, os, sys, numpy
from scipy.linalg import block_diag
from numpy.matlib import zeros

from math3d import Vector

from ..structures.atom_vector import AtomVector, AtomKeys
from ..structures.cell import Cell
from turbo_format import TurboTemplate
from turbo_basis import TurboBasis
from mos_format import MosReader


class ControlFormat:
	def __init__(self):
		self.base_format = None
		self.cell = None
		self.dm = []

	def species_map(self):
		ap = self.base_format.param('atoms')
		n = 0
		smap = {}
		while n < len(ap.multiparam):
			nums = set()
			ls = ap.multiparam[n].split()
			for r in ls[1].split(','):
				lls = r.split('-')
				if len(lls) == 1:
					nums.add(int(lls[0]))
				elif len(lls) == 2:
					nums |= set(range(int(lls[0]), int(lls[1]) + 1))
			n += 1
			name = ''
			while ap.multiparam[n - 1][-1] == '\\':
				l = re.split("(\s+\\\\)", ap.multiparam[n].strip())[0]
				ls = re.split("(\s+=)", l)
				if ls[0] == 'basis':
					name = ls[2]
				n += 1
			if name:
				for nn in nums:
					smap[nn] = name
		return smap

	def dft_matrix(self, n, na, nb):
		mat = zeros((n, n))
		for i in range(na):
			mat[i*2, i*2] = 1
		for i in range(nb):
			mat[i*2 + 1, i*2 + 1] = 1
		for i in range(n):
			print mat[i, i]
		return mat

	def load(self, of):
		self.base_format = of
		bp = self.base_format.param('basis')
		assert bp.lineparam.strip() == 'file=basis'
		self.bases = TurboBasis.read_basis('basis')
		smap = self.species_map()
		cp = self.base_format.param('coord')
		assert cp.lineparam.strip() == 'file=coord'
		coordf = TurboTemplate.from_file('coord')
		cp = coordf.param('coord')
		atoms = []
		for n, l in enumerate(cp.multiparam):
			ls = l.split()
			v = Vector(float(ls[0]), float(ls[1]), float(ls[2]))
			name = ls[3]
			a = AtomVector(name, v)
			atoms.append(a)
			a.data()[AtomKeys.ORBITAL_COUNT] = self.bases[smap[n + 1]].orbnum()
			a.data()[AtomKeys.ORBITAL_ARRAY] = self.bases[smap[n + 1]].orbarray()
		self.cell = Cell(atoms, [])

		twecp = self.base_format.param('twocomp-ecp')
		if twecp:
			ts = self.base_format.param('twocomp')
			ls = re.split("[\s\-]+", ts.multiparam[0].strip())
			nocc = len(range(int(ls[1]), int(ls[2]) + 1))

			rm = self.base_format.param('uhfmo_real')
			im = self.base_format.param('uhfmo_imag')
			assert rm.lineparam.strip() == 'file=realmos'
			assert im.lineparam.strip() == 'file=imagmos'
			rmat = MosReader.from_file('realmos').matrix
			imat = MosReader.from_file('imagmos').matrix
			basis_mat = []
			N = len(rmat) / 2
			for i in range(N * 2):
				tm = []
				for j in range(N):
					e = rmat[i][j] + imat[i][j] * 1j
					tm.append(e)
					e = rmat[i][j + N] + imat[i][j + N] * 1j
					tm.append(e)
				basis_mat.append(tm)
			bm = numpy.matrix(basis_mat)
#			occ_m = self.dft_matrix(N * 2, nalpha, nbeta)
			occ_m = block_diag(numpy.identity(nocc), numpy.zeros((N * 2 - nocc, N * 2 - nocc)))
			self.dm = bm.getH().dot(occ_m).dot(bm)
		else:
			assert False #TODO
			nalpha = int(re.split("[\s\-]+", self.base_format.param('alpha').multiparam[0].strip())[2])
			nbeta = int(re.split("[\s\-]+", self.base_format.param('beta').multiparam[0].strip())[2])
			no = self.base_format.param('natural')
			assert no.lineparam.split()[-1] == 'file=natorb'
			mat = MosReader.from_file('natorb').matrix
			basis_mat = []
			N = len(mat)
			for i in range(N * 2):
				tm = []
				for j in range(N * 2):
					e = 0
					if (i-j) % 2 == 0:
						e = mat[i /2 ][j / 2]
					tm.append(e)
				basis_mat.append(tm)
			bm = numpy.matrix(basis_mat)
			occ_m = self.dft_matrix(N * 2, nalpha, nbeta)
			self.dm = bm.getH().dot(occ_m).dot(bm)

	def get_format(self):
		return self.base_format

	def to_string(self):
		return self.get_format().to_string()

	def to_file(self, name):
		self.get_format().to_file(name)

	@staticmethod
	def from_string(datastring):
		res = ControlFormat()
		res.load(TurboTemplate.from_string(datastring))
		return res

	@staticmethod
	def from_file(name):
		res = ControlFormat()
		res.load(TurboTemplate.from_file(name))
		return res

