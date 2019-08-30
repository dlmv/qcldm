import os, logging
from ..structures.bond_system import MullikenOverlapBondData, LinearSystemChargeTransferBondData, DumbBondData
from ..structures.cluster_comparator import compare_clusters
from ..util.xyz_format import write_xyz
from ..util.fileutils import make_dir
from ..util.units import Units
from ..structures.atom_vector import AtomKeys

CONTROL_TEMPLATE = '''
$title
%%NAME%%
$symmetry c1
$coord    file=coord_ca
$atoms
%%ATOMS%%
$pople   AO
$basis    file=basis
$ecp    file=basis
$scfiterlimit       99
$scfconv        7
$thize     0.10000000E-04
$thime        5
$scfdamp   start=1.500  step=0.050  min=0.100
$scfdump
$scfintunit
 unit=30       size=0        file=twoint
$scfdiis   start=0.5
$scforbitalshift  automatic=0.1
$drvopt
   cartesian  on
   basis      off
   global     off
   hessian    on
   dipole     on
   nuclear polarizability
$interconversion  off
   qconv=1.d-10
   maxiter=25
$optimize
   internal   off
   cartesian  on
   global     off
   basis      off   logarithm
$coordinateupdate
   dqmax=0.3
   interpolate  on
   statistics    5
$forceupdate
   ahlrichs numgeo=0  mingeo=3 maxgeo=4 modus=<g|dq> dynamic fail=0.1
   threig=0.005  reseig=0.005  thrbig=3.0  scale=1.00  damping=0.0
$forceinit on
   diag=default
$lock off
$twocomp-ecp
$dft-section
$dft-functional pbe0
$dft-gridtype cvw-3
$ncpus   28
$scfmo    file=mos
$twocomp shells
 a    1-  %%ELECTRONS%%    (  1  )
$uhfmo_real       file=realmos
$uhfmo_imag       file=imagmos
$closed shells
$end'''[1:]

class AtomDataPart:
	def __init__(self, name, basis, ecp):
		self.name = name
		self.numbers = []
		self.basis = basis
		self.ecp = ecp
	
	def numstr(self):
		lastfirst = self.numbers[0]
		current = self.numbers[0]
		res = ''
		for x in self.numbers[1:] + [999]:
			if current == x-1:
				current = x
			elif current == lastfirst:
				res += "%d," % current
				current = x
				lastfirst = x
			else:
				res += "%d-%d," % (lastfirst, current)
				lastfirst = x
				current = x
		return res[:-1]
	
	def string(self):
		ls = []
		headstr = '{:<3}'.format(self.name.lower())
		headstr += self.numstr()
		headstr = '{:<79}\\'.format(headstr)
		ls.append(headstr)
		basisstr = '   basis =%s %s' % (self.name.lower(), self.basis) if self.basis else '   basis =none'
		if self.ecp:
			basisstr = '{:<79}\\'.format(basisstr)
			ls.append(basisstr)
			ecpstr = '   ecp   =%s %s' % (self.name.lower(), self.ecp)
			ls.append(ecpstr)
		else:
			ls.append(basisstr)
		return "\n".join(ls)
		
class AtomData:
	def __init__(self):
		self.atoms = {}
		self.names = []

	def add(self, name, number, full, cluster):
		replacement = name
		if not full and name != 'q':
			replacement, core = cluster.settings.embedding_map[name]
		if (name, full) in self.atoms.keys():
			self.atoms[(name, full)].numbers.append(number)
		else:
			basis = cluster.settings.basis_map.get((name, full))
			ecp = cluster.settings.ecp_map.get((name, full))
			atom = AtomDataPart(replacement, basis, ecp)
			atom.numbers.append(number)
			self.atoms[(name, full)] = atom
			self.names.append((name, full))
	
	def string(self):
		return "\n".join([self.atoms[x].string() for x in self.names])

class TurboWriter:

	@staticmethod
	def atoms_string(cluster):
		data = AtomData()
		for i, ca in enumerate(cluster.atoms[:len(cluster.core_atoms)]):
			data.add(ca.origin.name(), i+1, True, cluster)
		for j, ba in enumerate(cluster.atoms[len(cluster.core_atoms):len(cluster.core_atoms) + len(cluster.border_atoms)]):
			data.add(ba.origin.name(), i+j+2, False, cluster)
		for k, ea in enumerate(cluster.atoms[len(cluster.core_atoms) + len(cluster.border_atoms):len(cluster.core_atoms) + len(cluster.border_atoms) + len(cluster.electrostatic_atoms)]):
			data.add('q', i+j+k+3, False, cluster)
		return data.string()
		

	@staticmethod
	def write_control(cluster):
		with open(os.path.join(cluster.settings.name, 'control'), 'w') as ctrlf:
			data = CONTROL_TEMPLATE
			data = data.replace('%%ELECTRONS%%', str(cluster.total_electrons()))
			data = data.replace('%%NAME%%', str(cluster.settings.name))
			data = data.replace('%%ATOMS%%', TurboWriter.atoms_string(cluster))
			ctrlf.write(data)
			
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

	@staticmethod
	def write_coord_ca(cluster):
		k = Units.UNIT / Units.BOHR
		with open(os.path.join(cluster.settings.name, 'coord_ca'), 'w') as cf:
			cf.write('$coord\n')
			for ca in cluster.atoms[:len(cluster.core_atoms)]:
				cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ca.origin.position().x * k, ca.origin.position().y * k, ca.origin.position().z * k, ca.origin.name()))
			for ba in cluster.atoms[len(cluster.core_atoms):len(cluster.core_atoms) + len(cluster.border_atoms)]:
				replacement, core = cluster.settings.embedding_map[ba.origin.name()]
				cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ba.origin.position().x * k, ba.origin.position().y * k, ba.origin.position().z * k, replacement))
			for ea in cluster.atoms[len(cluster.core_atoms) + len(cluster.border_atoms):len(cluster.core_atoms) + len(cluster.border_atoms) + len(cluster.electrostatic_atoms)]:
				cf.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(ea.origin.position().x * k, ea.origin.position().y * k, ea.origin.position().z * k, 'q'))
			cf.write('$end')

	@staticmethod
	def write_mos(cluster):
		with open(os.path.join(cluster.settings.name, 'mos'), 'w') as mosf:
			mosf.write("$scfmo    expanded   format(4d20.14)\n$end")






					








