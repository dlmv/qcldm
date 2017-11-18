import os, logging
from bond_system import MullikenOverlapBondData, LinearSystemChargeTransferBondData
from ..util.xyz_format import write_xyz
from ..util.fileutils import make_dir
from ..util.units import Units
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


	def write_embedding(self, dirname='.'):
		k = Units.UNIT / Units.BOHR
		if os.path.exists('embedding.template'):
			embeds = {}
			with open('embedding.template') as f:
				for l in f.read().splitlines():
					a, na, core = l.split()
					embeds[a] = [na, int(core)]
			with open(os.path.join(dirname, 'embedding'), 'w') as ef, open(os.path.join(dirname, 'embedding.start'), 'w') as esf, open(os.path.join(dirname, 'coord'), 'w') as cf:
				cf.write('$coord\n')
				for ca in self.atoms[:len(self.core_atoms)]:
					cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ca.origin.position().x * k, ca.origin.position().y * k, ca.origin.position().z * k, ca.origin.name()))
				for ba in self.atoms[len(self.core_atoms):len(self.core_atoms) + len(self.border_atoms)]:
					name, core = embeds[ba.origin.name()]
					ef.write("{:3}  {:7.3f}\n".format(name, ba.valence + core))
					esf.write("{:3}  {:7.3f}  {:7.3f}  {:7.3f}\n".format(name, ba.valence + core, core, core + ba.origin.data()[AtomKeys.ESTIMATED_VALENCE]))
					cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ba.origin.position().x * k, ba.origin.position().y * k, ba.origin.position().z * k, 'zz'))
				for ea in self.atoms[len(self.core_atoms) + len(self.border_atoms):len(self.core_atoms) + len(self.border_atoms) + len(self.electrostatic_atoms)]:
					ef.write("{:3}  {:7.3f}\n".format('q', ea.charge + core))
					cmin, cmax = min(0, ea.origin.data()[AtomKeys.BADER_CHARGE]), max(0, ea.origin.data()[AtomKeys.BADER_CHARGE])
					esf.write("{:3}  {:7.3f}  {:7.3f}  {:7.3f}\n".format('q', ea.charge, cmin, cmax))
					cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ea.origin.position().x * k, ea.origin.position().y * k, ea.origin.position().z * k, 'zz'))
				cf.write('$end')
					

	def write_structure(self, dirname='.'):
		make_dir(dirname)
		xyz_cluster = [ca.origin for ca in self.atoms]
		write_xyz(xyz_cluster, os.path.join(dirname, "cluster_structure.xyz"))

	def write_charges(self, dirname='.'):
		make_dir(dirname)
		res = '===============================================================\n' +\
			'             Valence       Charge    Original valence   Bader\n' + \
			'===============================================================\n'
		fmt = "  {:3}       {:7.3f}       {:7.3f}       {:7.3f}       {:7.3f}\n"
		tc = 0.
		tv = 0.
		for ca in self.atoms:
			full = ca.origin.data()[AtomKeys.ESTIMATED_VALENCE]
			bcharge = ca.origin.data()[AtomKeys.BADER_CHARGE]
			res += fmt.format(ca.origin.name(), ca.valence, ca.charge, full, bcharge)
			tc += ca.charge - ca.valence
			tv += ca.valence
		res += '===============================================================\n'
		res += 'Total:      {:7.3f}       {:7.3f}\n'.format(tv, tc)
		with open(os.path.join(dirname, "cluster.charge"), 'w') as f:
			f.write(res)







