import re, os, sys, math
from math3d import Vector
from openmx_format import Openmx_format
from openmx_species import openmx_species
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys

class DAT_INPUT:

	def __init__(self):
		self.base_format = None
		self.data_path = '.'
		self.system_name = '?'
		self.stdout_level = 1
		self.fileout_level = 1
		self.species = {}
		self.cell = None

	def get_unit(self, unit_str):
		if unit_str.lower() == 'au':
			return Units.BOHR / Units.UNIT
		elif unit_str.lower() == 'ang':
			return Units.ANGSTROM / Units.UNIT
		else:
			raise RuntimeError("Unknown unit: %s" % unit_str)

	def cur_unit_str(self):
		if Units.BOHR == Units.UNIT:
			return "AU"
		elif Units.ANGSTROM == Units.UNIT:
			return "Ang"
		else:
			raise RuntimeError("Unknown unit: %f" % Units.UNIT)

	def load(self, of):
		self.base_format = of

		curdir = self.base_format.param('System.CurrrentDirectory')[0] or "./"
		assert curdir == "./", "Cannot work with changed current dir"

		self.data_path = self.base_format.param('DATA.PATH')[0] or "../DFT_DATA13"
		self.system_name = self.base_format.param('System.Name')[0]
		self.stdout_level = int(self.base_format.param('level.of.stdout')[0])
		self.fileout_level = int(self.base_format.param('level.of.fileout')[0])

		spe_num = int(self.base_format.param('Species.Number')[0])
		assert spe_num == len(self.base_format.multiparam('Definition.of.Atomic.Species')), "Incorrect species number!"
		self.species = {}
		for ls in self.base_format.multiparam('Definition.of.Atomic.Species'):
			s = openmx_species(ls, self.data_path)
			self.species[s.name] = s

		atom_num = int(self.base_format.param('Atoms.Number')[0])
		assert atom_num == len(self.base_format.multiparam('Atoms.SpeciesAndCoordinates')), "Incorrect atom number!"
		atom_unit = self.get_unit(self.base_format.param('Atoms.SpeciesAndCoordinates.Unit')[0])
		atoms = []
		for ls in self.base_format.multiparam('Atoms.SpeciesAndCoordinates'):
			v = Vector(float(ls[2]) * atom_unit, float(ls[3]) * atom_unit, float(ls[4]) * atom_unit)
			name = self.species[ls[1]].name
			a = AtomVector(name, v)
			a.data()[AtomKeys.ORBITAL_COUNT] = self.species[a.name()].orbnum()
			a.data()[AtomKeys.ORBITAL_ARRAY] = self.species[a.name()].orbarray()
			a.data()[AtomKeys.CUTOFF] = self.species[a.name()].real_r()
			if self.species[a.name()].basis:
				a.data()[AtomKeys.FULL_VALENCE] = self.species[a.name()].basis.eval
			atoms.append(a)

		vectors = []
		if self.base_format.multiparam('Atoms.UnitVectors'):
			cell_unit = self.get_unit(self.base_format.param('Atoms.UnitVectors.Unit')[0])
			for ls in self.base_format.multiparam('Atoms.UnitVectors'):
				v = Vector(float(ls[0]) * cell_unit, float(ls[1]) * cell_unit, float(ls[2]) * cell_unit)
				vectors.append(v)
			assert len(vectors) == 3, "Wrong vector count!"
		self.cell = Cell(atoms, vectors)


	def get_format(self):

		self.base_format.set_param('System.CurrrentDirectory', ["./"])

		self.base_format.set_param('DATA.PATH', [self.data_path])
		self.base_format.set_param('System.Name', [self.system_name])
		self.base_format.set_param('level.of.stdout', [str(self.stdout_level)])
		self.base_format.set_param('level.of.fileout', [str(self.fileout_level)])
		self.base_format.set_param('Species.Number', [str(len(self.species.keys()))])
		ls = []
		for name in sorted(self.species.keys()):
			ls.append(self.species[name].fulllist())
		self.base_format.set_multiparam("Definition.of.Atomic.Species", ls)

		self.base_format.set_param('Atoms.Number', [str(len(self.cell.atoms))])
		self.base_format.set_param('Atoms.SpeciesAndCoordinates.Unit', [self.cur_unit_str()])
		ls = []
		for i in range(len(self.cell.atoms)):
			a = self.cell.atoms[i]
			sp = self.species[a.name()]
			occ = (sp.pp.eval + sp.pp.z - sp.pp.etotal) / 2.
			ls.append([str(i+1), a.name(), "{:19.14f}".format(a.position()[0]), "{:19.14f}".format(a.position()[1]), "{:19.14f}".format(a.position()[2]), "{:19.14f}".format(occ), "{:19.14f}".format(occ)])
		self.base_format.set_multiparam("Atoms.SpeciesAndCoordinates", ls)

		if self.cell.vectors:
			ls = []
			self.base_format.set_param('Atoms.UnitVectors.Unit', [self.cur_unit_str()])
			for v in self.cell.vectors:
				ls.append(["{:19.14f}".format(v[0]), "{:19.14f}".format(v[1]), "{:19.14f}".format(v[2])])
			self.base_format.set_multiparam("Atoms.UnitVectors", ls)
			self.base_format.set_param('scf.EigenvalueSolver', ['band'])
		else:
			self.base_format.remove_params(['Atoms.UnitVectors.Unit', "Atoms.UnitVectors"])
			self.base_format.set_param('scf.EigenvalueSolver', ['cluster'])
		return self.base_format

	def to_string(self):
		return self.get_format().to_string()

	def to_file(self, name):
		self.get_format().to_file(name)

	@staticmethod
	def from_string(datastring):
		res = DAT_INPUT()
		res.load(Openmx_format.from_string(datastring))
		return res

	@staticmethod
	def from_file(name):
		res = DAT_INPUT()
		res.load(Openmx_format.from_file(name))
		return res


