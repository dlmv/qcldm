import os, logging
from bond_system import MullikenOverlapBondData, LinearSystemChargeTransferBondData
from ..util.xyz_format import write_xyz
from ..util.fileutils import make_dir
from ..structures.atom_vector import AtomKeys


class ClusterAtom:

	def __init__(self, origin, valence, charge):
		self.origin = origin
		self.valence = valence
		self.charge = charge

		self.atoms = []

class Cluster:
	
	def __init__(self, cell, center, shells_num, electro_shells_num=0):
		self.cell = cell
		self.center = center
		self.shells_num = shells_num

		shells = cell.neighbours.neighbours_cluster(center, shells_num)

		border_shell_pos = shells_num - electro_shells_num
		assert border_shell_pos > 1, "Too few shells for a cluster!"

		self.core_atoms = []
		self.border_atoms = []
		self.electrostatic_atoms = []

		for i in range(len(shells)):
			for a in shells[i]:
				if i < border_shell_pos:
					self.core_atoms.append(a)
				elif i == border_shell_pos:
					self.border_atoms.append(a)
				else:
					self.electrostatic_atoms.append(a)

	def round_valence(self, atoms):
		logging.debug(u'Rounding total valence...')
		tv = 0
		for a in atoms:
			tv += a.valence
		itv = round(tv)
		for a in atoms:
			nv = a.valence * itv/tv
			a.charge += nv - a.valence
			a.valence = nv

	def rearrange_charges(self, atoms):
		logging.debug(u'Checking for exceeding values...')
		exc_atoms = set()
		exc_chg = 0.
		for a in atoms:
			valence = a.origin.data()[AtomKeys.ESTIMATED_VALENCE]
			if a.charge > valence:
				exc_chg += a.charge - valence
				exc_atoms.add(a.origin)
				a.charge = valence
		n_good = len(atoms) - len(exc_atoms)
		assert n_good > 0, "All atoms exceeding valence & charges!!!"
		ok = True
		for a in atoms:
			valence = a.origin.data()[AtomKeys.ESTIMATED_VALENCE]
			ok = ok and a.charge <= valence
			if a.origin not in exc_atoms:
				a.charge += exc_chg / n_good
		if not ok:
			self.rearrange_charges(atoms)


	def estimate_charges(self, dm, olp):
		self.ct_data = LinearSystemChargeTransferBondData(self.cell) if self.electrostatic_atoms else None
		self.mul_data = MullikenOverlapBondData(dm, olp)

		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Estimating cluster charges')
		logging.info(u'*********************************************')
		logging.info(u'')
		atoms = []
		logging.debug(u'Core atoms')
		tmp = []
		for i in range(len(self.core_atoms)):
			a = self.core_atoms[i]
			valence = a.data()[AtomKeys.ESTIMATED_VALENCE]
			ca = ClusterAtom(a, valence, valence)
			tmp.append(ca)
			logging.debug(u'  %d/%d' % (i + 1, 	len(self.core_atoms)))

		atoms.extend(tmp)
		tmp = []

		logging.debug(u'Border atoms')
		for i in range(len(self.border_atoms)):
			a = self.border_atoms[i]
			nbs = self.cell.neighbours.first_neighbours(a)
			total_overlap = 0
			inside_overlap = 0
			el_charge = 0
			for nb in nbs:
				molp = self.mul_data(a, nb)
				total_overlap += molp
				if nb in self.core_atoms or nb in self.border_atoms:
					inside_overlap += molp
				elif nb in self.electrostatic_atoms:
					add_charge = self.ct_data(a, nb)
					el_charge += add_charge
			valence = a.data()[AtomKeys.ESTIMATED_VALENCE]
			valence *= inside_overlap / total_overlap
			ca = ClusterAtom(a, valence, el_charge + valence)
			tmp.append(ca)
			logging.debug(u'  %d/%d' % (i + 1, len(self.border_atoms)))


		self.round_valence(tmp)
		self.rearrange_charges(tmp)
		atoms.extend(tmp)

		tmp = []
		logging.debug(u'Electrostatic atoms')
		for i in range(len(self.electrostatic_atoms)):
			a = self.electrostatic_atoms[i]
			nbs = self.cell.neighbours.first_neighbours(a)
			el_charge = 0
			for nb in nbs:
				if nb in self.core_atoms or nb in self.border_atoms or nb in self.electrostatic_atoms:
					add_charge = self.ct_data(a, nb)
					el_charge += add_charge
			ca = ClusterAtom(a, 0, el_charge)
			tmp.append(ca)
			logging.debug(u'  %d/%d' % (i + 1, len(self.electrostatic_atoms)))

		atoms.extend(tmp)

		self.atoms = atoms

	def write_structure(self, dirname='.'):
		make_dir(dirname)
		xyz_cluster = [ca.origin for ca in self.atoms]
		write_xyz(xyz_cluster, os.path.join(dirname, "cluster.xyz"))

	def write_charges(self, dirname='.'):
		make_dir(dirname)
		res = '=========================================================\n' +\
			'             Valence      Charge      Original valence\n' + \
			'=========================================================\n'
		fmt = "  {:3}       {:7.3f}       {:7.3f}       {:7.3f}\n"
		tc = 0.
		tv = 0.
		for ca in self.atoms:
			full = ca.origin.data()[AtomKeys.ESTIMATED_VALENCE]
			res += fmt.format(ca.origin.name(), ca.valence, ca.charge, full)
			tc += ca.charge - ca.valence
			tv += ca.valence
		res += '=========================================================\n'
		res += 'Total:      {:7.3f}       {:7.3f}\n'.format(tv, tc)
		with open(os.path.join(dirname, "cluster.charge"), 'w') as f:
			f.write(res)







