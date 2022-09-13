import re, os, sys, math
from .openmx_format import Openmx_format
from ..atom.pseudo_potential import SeparablePseudoPotential
from ..functions.numeric_function import NumericFunction


class VPS:

	@staticmethod
	def pp_from_format(of):
		pp = SeparablePseudoPotential()
		pp.z = float(of.param('AtomSpecies')[0])
		pp.etotal = float(of.param('total.electron')[0])
		pp.eval = float(of.param('valence.electron')[0])

		jON = of.param('j.dependent.pseudo.potentials')[0].lower() == 'on'
		assert jON, "non-relativistic not implemented yet(?)"

		pccON = of.param('charge.pcc.calc')[0].lower() == 'on'

		projN = int(of.multiparam('project.energies')[0][0])
		assert len(of.multiparam('project.energies')[1:]) == projN
		proj_energies = []
		proj_ls = []
		projectors = []
		for ls in of.multiparam('project.energies')[1:]:
			proj_ls.append(int(ls[0]))
			proj_energies.append((float(ls[1]), float(ls[2])))

		pp.loc = NumericFunction()
		for i in range(projN):
			projectors.append((NumericFunction(), NumericFunction()))
		for ls in of.multiparam('Pseudo.Potentials'):
			assert len(ls) == projN * 2 + 3
			r = float(ls[1])
			pp.loc.data.append((r, float(ls[2])))
			for i in range(projN):
				projectors[i][0].data.append((r, float(ls[i * 2 + 3])))
				projectors[i][1].data.append((r, float(ls[i * 2 + 4])))
		pp.blochl_projectors = list(zip(proj_ls, projectors, proj_energies))

		if pccON:
			pp.pcc = NumericFunction()
			for ls in of.multiparam('density.PCC'):
				r = float(ls[1])
				pp.pcc.data.append((r, float(ls[2])))

		pp.grid = pp.loc.grid()
		return pp


	@staticmethod
	def format_from_pp(pp):
		res = Openmx_format([])
		res.set_param('AtomSpecies', ["{}".format(pp.z)])
		res.set_param('total.electron', ["{}".format(pp.etotal)])
		res.set_param('valence.electron', ["{}".format(pp.eval)])
		res.add_text('')
		res.set_param('grid.num.output', ["{}".format(len(pp.grid))])
		res.add_text('')
		res.set_param('charge.pcc.calc', ["on" if pp.pcc else "off"])
		res.add_text('')
		res.set_param('j.dependent.pseudo.potentials', ["on"])
		res.add_text('')
		ls = [[str(len(pp.blochl_projectors))]]
		for l, (p1, p2), (e1, e2) in pp.blochl_projectors:
			tls = [str(l), "{:20.14e}".format(e1), "{:20.14e}".format(e2)]
			ls.append(tls)
		res.set_multiparam("project.energies", ls)
		res.add_text('')
		ls = []
		for i, r in enumerate(pp.grid):
			tls = ["{:20.14e}".format(math.log(r)), "{:20.14e}".format(r), "{:20.14e}".format(pp.loc.data[i][1])]
			for l, (p1, p2), (e1, e2) in pp.blochl_projectors:
				tls.append("{:20.14e}".format(p1.data[i][1]))
				tls.append("{:20.14e}".format(p2.data[i][1]))
			ls.append(tls)
		res.set_multiparam("Pseudo.Potentials", ls)
		res.add_text('')
		ls = []
		if pp.pcc:
			for r, p in pp.pcc.data:
				tls = ["{:20.14e}".format(math.log(r)), "{:20.14e}".format(r), "{:20.14e}".format(p)]
				ls.append(tls)
			res.set_multiparam("density.PCC", ls)
		return res

	@staticmethod
	def to_string(pp):
		return VPS.format_from_pp(pp).to_string()

	@staticmethod
	def to_file(pp, name):
		return VPS.format_from_pp(pp).to_file(name)

	@staticmethod
	def from_string(datastring):
		return VPS.pp_from_format(Openmx_format.from_string(datastring))

	@staticmethod
	def from_file(name):
		return VPS.pp_from_format(Openmx_format.from_file(name))


