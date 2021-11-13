import re, os, sys, math
from openmx_format import Openmx_format
from ..atom.basis import Basis
from ..functions.numeric_function import NumericFunction


class PAO:

	@staticmethod
	def basis_from_format(of):
		basis = Basis()
		basis.z = float(of.param('AtomSpecies')[0])
		basis.etotal = float(of.param('total.electron')[0])
		basis.eval = float(of.param('valence.electron')[0])
		basis.cutoff = float(of.param('radial.cutoff.pao')[0])
		lMax = int(of.param('PAO.Lmax')[0])
		mul = int(of.param('PAO.Mul')[0])

		basis.vden = NumericFunction()
		for ls in of.multiparam('valence.charge.density'):
			r = float(ls[1])
			basis.vden.data.append((r, float(ls[2])))

		basis.functions = {}
		for l in range(lMax + 1):
			lpaos = {}
			for i in range(mul):
				f = NumericFunction()
				f.l = l
				f.n = i + 1
				lpaos[f.n] = f
			for ls in of.multiparam('pseudo.atomic.orbitals.L={}'.format(l)):
				assert len(ls) == mul + 2
				r = float(ls[1])
				for i in range(mul):
					lpaos[i + 1].data.append((r, float(ls[i + 2])))
			basis.functions[l] = lpaos

		basis.grid = basis.vden.grid()
		return basis

	@staticmethod
	def format_from_basis(basis):
		basis.align_functions()
		res = Openmx_format([])
		res.set_param('AtomSpecies', ["{}".format(basis.z)])
		res.set_param('total.electron', ["{}".format(basis.etotal)])
		res.set_param('valence.electron', ["{}".format(basis.eval)])
		res.add_text('')
		res.set_param('grid.num.output', ["{}".format(len(basis.grid))])
		res.add_text('')
		res.set_param('radial.cutoff.pao', ["{}".format(basis.cutoff)])
		res.add_text('')
		ls = []
		if not basis.vden:
			basis.vden = basis.estimate_full_vden()
		for r, v in basis.vden.data:
			tls = ["{:18.15f}".format(math.log(r)), "{:18.15f}".format(r), "{:18.15f}".format(v)]
			ls.append(tls)
		res.set_multiparam("valence.charge.density", ls)
		res.add_text('')
		res.set_param('PAO.Lmax', ["{}".format(len(basis.functions.keys()) - 1)])
		res.set_param('PAO.Mul', ["{}".format(len(basis.functions.values()[0]))])
		for l in sorted(basis.functions.keys()):
			lpaos = basis.functions[l]
			ls = []
			for i, r in enumerate(basis.grid):
				tls = ["{:18.15f}".format(math.log(r)), "{:18.15f}".format(r)]
				for n in sorted(lpaos.keys()):
					lp = lpaos[n]
					tls.append("{:18.15f}".format(lp.data[i][1]))
				ls.append(tls)
			res.set_multiparam('pseudo.atomic.orbitals.L={}'.format(l), ls)

		return res

	@staticmethod
	def to_string(basis):
		return PAO.format_from_basis(basis).to_string()

	@staticmethod
	def to_file(basis, name):
		return PAO.format_from_basis(basis).to_file(name)

	@staticmethod
	def from_string(datastring):
		return PAO.basis_from_format(Openmx_format.from_string(datastring))

	@staticmethod
	def from_file(name):
		return PAO.basis_from_format(Openmx_format.from_file(name))
