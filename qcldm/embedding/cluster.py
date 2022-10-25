import os, logging
from ..structures.bond_system import MullikenOverlapBondData, LinearSystemChargeTransferBondData, DumbBondData
from ..structures.cluster_comparator import compare_clusters
from ..util.xyz_format import write_xyz
from ..util.fileutils import make_dir
from ..util.units import Units
from ..structures.atom_vector import AtomKeys
from ..atom.shells import Shells
from .turbo_writer import TurboWriter


class ClusterAtom:

	def __init__(self, origin, valence, charge):
		self.origin = origin
		self.valence = valence
		self.charge = charge

		self.atoms = []



class Cluster:

	def get_add_center(self, centers, skip_centers, name, shell, r):
		center = centers[0]
		shells = self.cell.neighbours.neighbours_cluster([center], shell, self.settings.bond_distance_override_map)
		for a in shells[-1]:
			if a in centers or a in skip_centers:
				continue
			if a.name() == name and abs(center.distance(a) - r) < 1e-4:
				return a
		assert False, 'additional center not found! {} {} {}' .format(name, shell, r)

	def add_cations(self, shells, name, r):
		for a in self.cell.extended_cell(2):
			if a.name() == name and a.distance(self.centers[0]) <= r:
				print(a)
				found = False
				for shell in shells:
					if a in shell:
						found = True
				if not found:
					shells[-1].append(a)

	def __init__(self, cell, settings):
		self.cell = cell
		self.settings = settings
		
		center = cell.cell[settings.center - 1]
		assert center.name().lower() == settings.centername.lower(), 'wrong center name: %s %s' % (center.name(), settings.centername)
		self.centers = [center]
		
		skip_centers = []
		
		for (atom, shell, r, add_not_skip) in self.settings.add_centers:
			ac = self.get_add_center(self.centers, skip_centers, atom, shell, r)
			if add_not_skip:
				self.centers.append(ac)
			else:
				skip_centers.append(ac)
		
		self.estimate_atoms_charges()

		shells = cell.neighbours.neighbours_cluster(self.centers, self.settings.inner_shell_num, self.settings.bond_distance_override_map)
		
		if self.settings.expand_covalent:
			cell.neighbours.expand_covalent_border(shells, self.settings.bond_distance_override_map)
		
		cell.neighbours.expand_neighbours(shells, self.settings.bond_distance_override_map)
		
		for name, r in self.settings.cation_sphere:
			self.add_cations(shells, name, r)
		
		for i in range(settings.electro_shell_num):
			cell.neighbours.expand_neighbours(shells, self.settings.bond_distance_override_map)

#		assert settings.inner_shell_num > 0, "Too few shells for a cluster!"

		self.core_atoms = []
		self.border_atoms = []
		self.electrostatic_atoms = []

		for i in range(len(shells)):
			for a in shells[i]:
				if i < settings.inner_shell_num + 1:
					self.core_atoms.append(a)
				elif i == settings.inner_shell_num + 1:
					self.border_atoms.append(a)
				else:
					self.electrostatic_atoms.append(a)
#		self.core_atoms.sort()
#		self.border_atoms.sort()
#		self.electrostatic_atoms.sort()

	def estimate_atoms_charges(self):
		for a in self.cell.atoms:
			if not AtomKeys.ESTIMATED_VALENCE in list(a.data().keys()):
				a.data()[AtomKeys.ESTIMATED_VALENCE] = Shells.estimate_valence_byname(a.name())
		
			if not AtomKeys.ESTIMATED_CHARGE in list(a.data().keys()):
				a.data()[AtomKeys.ESTIMATED_CHARGE] = Shells.estimate_charge_byname(a.name())
			
			if a.name() in list(self.settings.valence_override_map.keys()):
				a.data()[AtomKeys.ESTIMATED_VALENCE] = self.settings.valence_override_map[a.name()]
			if a.name() in list(self.settings.charge_override_map.keys()):
				a.data()[AtomKeys.ESTIMATED_CHARGE] = self.settings.charge_override_map[a.name()]
		

	def round_valence(self, atoms, desired):
		logging.debug('Rounding total valence...')
		tv = 0
		for a in atoms:
			tv += a.valence
		logging.debug('  Current is %f; rounding to %d' % (tv, desired))
#		itv = round(tv)
		itv = desired
		for a in atoms:
			nv = a.valence * itv/tv
			a.charge += nv - a.valence
			a.valence = nv

	def estimate_central_charge(self):
		s = 0
		for a in self.core_atoms:
			s += a.data()[AtomKeys.ESTIMATED_CHARGE]
		return s

	def total_electrons(self):
		s = -self.estimate_central_charge()
		for a in self.core_atoms:
			s += a.data()[AtomKeys.FULL_VALENCE]
		return s
		

	def rearrange_charges(self, atoms):
		logging.debug('Checking for exceeding values...')
		add_caps = []
		sub_caps = []
		exc_chg = 0.
		for a in atoms:
			valence = a.origin.data()[AtomKeys.ESTIMATED_VALENCE]
			if a.charge > valence:
				exc_chg += a.charge - valence
				a.charge = valence
				add_caps.append(0)
				sub_caps.append(valence)
			elif a.charge < 0:
				exc_chg += a.charge
				a.charge = 0
				add_caps.append(valence)
				sub_caps.append(0)
			else:
				add_caps.append(valence - a.charge)
				sub_caps.append(a.charge)
		if exc_chg > 0:
			s = sum(add_caps)
			logging.debug('  Dustributing positive charge of %f' % exc_chg)
			if s < exc_chg:
				logging.error('  But only %f can be distributed!' % s)
				assert False
			for a, c in zip(atoms, add_caps):
				a.charge += exc_chg * c / s
		elif exc_chg < 0:
			exc_chg *= -1
			s = sum(sub_caps)
			logging.debug('  Dustributing negative charge of %f' % exc_chg)
			if s < exc_chg:
				logging.error('  But only %f can be distributed!' % s)
				assert False
			for a, c in zip(atoms, sub_caps):
				a.charge -= exc_chg * c / s
		else:
			logging.debug('  Everything already ok')


	def estimate_charges_mulliken(self, dm, olp,):
		self.estimate_charges(MullikenOverlapBondData(dm, olp))

	def estimate_charges_dumb(self):
		self.estimate_charges(DumbBondData())

	def estimate_charges(self, bond_data):
		self.charge_key = AtomKeys.MULLIKEN_CHARGE
		if self.settings.charge_method == 'bader':
			self.charge_key = AtomKeys.BADER_CHARGE
		elif self.settings.charge_method == 'redox':
			self.charge_key = AtomKeys.ESTIMATED_CHARGE
		self.ct_data = LinearSystemChargeTransferBondData(self.cell, self.charge_key, self.settings.bond_distance_override_map) if self.electrostatic_atoms else None
		self.mul_data = bond_data

		logging.info('')
		logging.info('*********************************************')
		logging.info('  Estimating cluster charges')
		logging.info('*********************************************')
		logging.info('')
		atoms = []
		logging.debug('Core atoms')
		tmp = []
		for i in range(len(self.core_atoms)):
			a = self.core_atoms[i]
			valence = a.data()[AtomKeys.ESTIMATED_VALENCE]
			ca = ClusterAtom(a, valence, valence)
			tmp.append(ca)
			logging.debug('  %d/%d' % (i + 1, 	len(self.core_atoms)))

		atoms.extend(tmp)
		desired = -self.estimate_central_charge()
		tmp = []

		logging.debug('Border atoms')
		for i in range(len(self.border_atoms)):
			a = self.border_atoms[i]
			nbs = self.cell.neighbours.first_neighbours(a, self.settings.bond_distance_override_map)
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
			logging.debug('  %d/%d' % (i + 1, len(self.border_atoms)))


		self.round_valence(tmp, desired)
		self.rearrange_charges(tmp)
		atoms.extend(tmp)

		tmp = []
		logging.debug('Electrostatic atoms')
		for i in range(len(self.electrostatic_atoms)):
			a = self.electrostatic_atoms[i]
			nbs = self.cell.neighbours.first_neighbours(a, self.settings.bond_distance_override_map)
			el_charge = 0
			for nb in nbs:
				if nb in self.core_atoms or nb in self.border_atoms or nb in self.electrostatic_atoms:
					add_charge = self.ct_data(a, nb)
					el_charge += add_charge
			ca = ClusterAtom(a, 0, el_charge)
			tmp.append(ca)
			logging.debug('  %d/%d' % (i + 1, len(self.electrostatic_atoms)))

		atoms.extend(tmp)

		self.atoms = atoms
		
	def make_groups(self):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  GROUPING ATOMS')
		logging.info('*********************************************')
		logging.info('')
		groups = {}
		atoms = self.border_atoms + self.electrostatic_atoms
		full_atoms = self.core_atoms + self.border_atoms + self.electrostatic_atoms
		index = 0
		for i in range(len(atoms)):
			logging.debug('  %d/%d' % (i + 1, len(atoms)))
			a = atoms[i]
			if a.tuple_data() in list(groups.keys()):
				continue
			found = False
			if self.settings.use_symmetry:	
				for j in range(len(atoms)):
					a1 = atoms[j]
					if a1.tuple_data() in list(groups.keys()) or i==j:
						continue
					if compare_clusters(a, full_atoms, a1, full_atoms):
						found = True
						groups[a.tuple_data()] = "GROUP_%d" % (index+1)
						groups[a1.tuple_data()] = "GROUP_%d" % (index+1)
			if not found:
				groups[a.tuple_data()] = "SINGLE_%d" % (index+1)
			index += 1
		return groups
							
	def write(self):
		self.write_structure()
		self.write_charges()
		if self.settings.make_turbo:
			self.write_embedding()
		
			
		


	def write_embedding(self):
		TurboWriter.write_embedding(self)
		TurboWriter.write_embedding_start(self)
		TurboWriter.write_coord(self)
		TurboWriter.write_coord_ca(self)
		TurboWriter.write_mos(self)
		TurboWriter.write_control(self)
					

	def write_structure(self):
		make_dir(self.settings.name)
		xyz_cluster = [ca.origin for ca in self.atoms]
		write_xyz(xyz_cluster, os.path.join(self.settings.name, "cluster_structure.xyz"))

	def write_charges(self):
		make_dir(self.settings.name)
		res = '===============================================================\n' +\
			'             Valence       Charge    Original valence   %s\n' % self.charge_key + \
			'===============================================================\n'
		fmt = "  {:3}       {:7.3f}       {:7.3f}       {:7.3f}       {:7.3f}\n"
		tc = 0.
		tv = 0.
		for ca in self.atoms:
			full = ca.origin.data()[AtomKeys.ESTIMATED_VALENCE]
			bcharge = ca.origin.data()[self.charge_key]
			res += fmt.format(ca.origin.name(), ca.valence, ca.charge, full, bcharge)
			tc += ca.charge - ca.valence
			tv += ca.valence
		res += '===============================================================\n'
		res += 'Total:      {:7.3f}       {:7.3f}\n'.format(tv, tc)
		with open(os.path.join(self.settings.name, "cluster.charge"), 'w') as f:
			f.write(res)







