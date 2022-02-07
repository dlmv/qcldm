import re, os, sys, numpy
from scipy.linalg import block_diag
from numpy.matlib import zeros

from math3d import Vector

from ..structures.atom_vector import AtomVector, AtomKeys
from ..structures.cell import Cell
from turbo_format import TurboTemplate
from turbo_basis import TurboBasis
from mos_format import MosReader

from ..util.elements import ELEMENTS
from ..util.units import Units


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
			basis_name = None
			ecp_name = None
			while ap.multiparam[n - 1][-1] == '\\':
				l = re.split("(\s+\\\\)", ap.multiparam[n].strip())[0]
				ls = re.split("(\s+=)", l)
				if ls[0] == 'basis':
					basis_name = ls[2]
				elif ls[0] == 'ecp':
					ecp_name = ls[2]
				n += 1
			for nn in nums:
				smap[nn] = (basis_name, ecp_name)
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

	@staticmethod
	def read_embedding():
		res = []
		if os.path.exists('embedding'):
			with open('embedding') as em:
				for line in em.read().splitlines():
					ls = line.split()
					res.append((ls[0], float(ls[1])))
		return res
					
		
	@staticmethod
	def get_normal_name(a):
		return a[0].upper() + a[1:].lower()

	def load(self, of, embmap):
		self.base_format = of
		bp = self.base_format.param('basis')
		assert bp.lineparam.strip() == 'file=basis'
		self.bases, self.ecps = TurboBasis.read_basis('basis')
		smap = self.species_map()
		
		embedding = ControlFormat.read_embedding()
			
		cp = self.base_format.param('coord')
		assert cp.lineparam.strip() == 'file=coord'
		coordf = TurboTemplate.from_file('coord')
		cp = coordf.param('coord')
		atoms = []
		empos = 0
		for n, l in enumerate(cp.multiparam):
			ls = l.split()
			k =  Units.BOHR / Units.UNIT
			v = Vector(float(ls[0]) * k, float(ls[1]) * k, float(ls[2]) * k)
			name = ControlFormat.get_normal_name(ls[3])
			
			embedded = False
			emcharge = 0
			
			if name.lower() == 'zz':
				name = embedding[empos][0]
				emcharge = embedding[empos][1]
				embedded = True
				empos += 1
			
			name = name[0].upper() + name[1:].lower()
			
			if name == 'Q':
				nelec = 0
			else:
				nelec = ELEMENTS[name].number
			
			a = AtomVector(name, v)
			
			if embedded:
				if name in embmap.keys():
					name = embmap[name]
					a = AtomVector(name, v)
				ncore = self.ecps[smap[n + 1][1]].ncore if smap[n + 1][1] else 0
#				a.data()[AtomKeys.ESTIMATED_VALENCE] = emcharge - self.ecps[smap[n + 1][1]].ncore
				a.data()[AtomKeys.ESTIMATED_CHARGE] = emcharge - ncore
				a.data()[AtomKeys.FULL_VALENCE] = 0
			else:
				a.data()[AtomKeys.FULL_VALENCE] = nelec - self.ecps[smap[n + 1][1]].ncore if smap[n + 1][1] else nelec
				
			
			a.data()[AtomKeys.ORBITAL_COUNT] = self.bases[smap[n + 1][0]].orbnum() if smap[n + 1][0] != 'none' else 0
			a.data()[AtomKeys.ORBITAL_ARRAY] = self.bases[smap[n + 1][0]].orbarray() if smap[n + 1][0] != 'none' else []
			
			
			atoms.append(a)

			

				
		self.cell = Cell(atoms, [], [], [], [])
		
	def unused(self):

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
	def from_string(datastring, embmap = {}):
		res = ControlFormat()
		res.load(TurboTemplate.from_string(datastring), embmap)
		return res

	@staticmethod
	def from_file(name, embmap = {}):
		res = ControlFormat()
		res.load(TurboTemplate.from_file(name), embmap)
		return res

