import os, logging
from ..structures.bond_system import MullikenOverlapBondData, LinearSystemChargeTransferBondData, DumbBondData
from ..structures.cluster_comparator import compare_clusters
from ..util.xyz_format import write_xyz
from ..util.fileutils import make_dir
from ..util.units import Units
from ..structures.atom_vector import AtomKeys

class TurboWriter:

	@staticmethod
	def write_embedding(cluster):
		with open(os.path.join(cluster.settings.name, 'embedding'), 'w') as ef:
		
			for ba in cluster.atoms[len(cluster.core_atoms):len(cluster.core_atoms) + len(cluster.border_atoms)]:
				replacement, core = cluster.settings.embedding_map[ba.origin.name()]
				ef.write("{:3}  {:9.5f}\n".format(replacement, ba.charge + core))
				
			for ea in cluster.atoms[len(cluster.core_atoms) + len(cluster.border_atoms):len(cluster.core_atoms) + len(cluster.border_atoms) + len(cluster.electrostatic_atoms)]:
				ef.write("{:3}  {:9.5f}\n".format('q', ea.charge))
				
	@staticmethod
	def write_embedding_start(cluster):
		groups = cluster.make_groups()
		with open(os.path.join(cluster.settings.name, 'embedding.start'), 'w') as esf:
		
			for ba in cluster.atoms[len(cluster.core_atoms):len(cluster.core_atoms) + len(cluster.border_atoms)]:
				replacement, core = cluster.settings.embedding_map[ba.origin.name()]
				g = '' if ba.origin.tuple_data() not in groups.keys() else groups[ba.origin.tuple_data()]
				esf.write("{:3}  {:9.5f}  {:9.5f}  {:9.5f} {}\n".format(replacement, ba.charge + core, core, core + ba.origin.data()[AtomKeys.ESTIMATED_VALENCE], g))
				
			for ea in cluster.atoms[len(cluster.core_atoms) + len(cluster.border_atoms):len(cluster.core_atoms) + len(cluster.border_atoms) + len(cluster.electrostatic_atoms)]:
				g = '' if ea.origin.tuple_data() not in groups.keys() else groups[ea.origin.tuple_data()]
				cmin, cmax = min(0, ea.origin.data()[cluster.charge_key]), max(0, ea.origin.data()[cluster.charge_key])
				esf.write("{:3}  {:9.5f}  {:9.5f}  {:9.5f} {}\n".format('q', ea.charge, cmin, cmax, g))
				
	@staticmethod
	def write_coord(cluster):
		k = Units.UNIT / Units.BOHR
		with open(os.path.join(cluster.settings.name, 'coord'), 'w') as cf:
			cf.write('$coord\n')
			for ca in cluster.atoms[:len(cluster.core_atoms)]:
				cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ca.origin.position().x * k, ca.origin.position().y * k, ca.origin.position().z * k, ca.origin.name()))
			for bea in cluster.atoms[len(cluster.core_atoms):len(cluster.core_atoms) + len(cluster.border_atoms) + len(cluster.electrostatic_atoms)]:
				cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(bea.origin.position().x * k, bea.origin.position().y * k, bea.origin.position().z * k, 'zz'))
			cf.write('$end')


					








