import os, logging
class EmbeddingSettings:
	def __init__(self):
		self.center = -1
		self.centername = ''
		self.inner_shell_num = 0
		self.electro_shell_num = 0
		self.expand_covalent = False
		self.embedding_map = {}
		self.charge_override_map = {}
		self.valence_override_map = {}
		self.bond_distance_override_map = {}
		self.basis_map = {}
		self.ecp_map = {}
		self.name = ''
		self.make_turbo = False
		self.charge_method = 'mulliken'
		self.add_centers = []
		self.cation_sphere = []
		self.reverse_embedding_map = {}
		self.ignore_symmetry = False

	@staticmethod
	def from_file(name):
		with open(name) as f:
			es = EmbeddingSettings()
			filename, file_extension = os.path.splitext(name)
			assert filename
			assert file_extension
			es.name = filename
			for line in f.read().splitlines():
				line = line.split('#')[0]
				if line:
					ls = line.strip().split()
					param = ls[0]
					if param == 'center':
						centernum, centername = read_2params(ls)
						es.center = int(centernum)
						es.centername = centername
					elif param == 'inner_shells':
						es.inner_shell_num = int(read_param(ls))
					elif param == 'electro_shells':
						es.electro_shell_num = int(read_param(ls))
					elif param == 'embedding':
						atom, replacement, core = read_3params(ls)
						core = int(core)
						es.embedding_map[atom] = [replacement, core]
					elif param == 'make_turbo':
						es.make_turbo = True
					elif param == 'expand_covalent':
						es.expand_covalent = True
					elif param == 'basis_full':
						atom, basis = read_2params(ls)
						es.basis_map[(atom, True)] = basis
					elif param == 'basis_border':
						atom, basis = read_2params(ls)
						es.basis_map[(atom, False)] = basis
					elif param == 'ecp_full':
						atom, ecp = read_2params(ls)
						es.ecp_map[(atom, True)] = ecp
					elif param == 'ecp_border':
						atom, ecp = read_2params(ls)
						es.ecp_map[(atom, False)] = ecp
					elif param == 'charge_method':
						es.charge_method = read_param(ls)
					elif param == 'override_bond':
						bond, length = read_2params(ls)
						length = float(length)
						bond = tuple(sorted(bond.split('-')))
						es.bond_distance_override_map[bond] = length
					elif param == 'override_charge':
						atom, value = read_2params(ls)
						value = int(value)
						es.charge_override_map[atom] = value
					elif param == 'override_valence':
						atom, value = read_2params(ls)
						value = int(value)
						es.valence_override_map[atom] = value
					elif param == 'skip_center':
						atom, shell, r = read_3params(ls)
						shell = int(shell)
						r = float(r)
						es.add_centers.append((atom, shell, r, False))
					elif param == 'add_center':
						atom, shell, r = read_3params(ls)
						shell = int(shell)
						r = float(r)
						es.add_centers.append((atom, shell, r, True))
					elif param == 'add_cations_in_sphere':
						atom, r = read_2params(ls)
						r = float(r)
						es.cation_sphere.append((atom, r))
					elif param == 'reverse_embedding':
						atom, replacement = read_2params(ls)
						core = int(core)
						es.reverse_embedding_map[atom] = replacement
					elif param == 'ignore_symmetry':
						es.ignore_symmetry = True
					else:
						logging.error('Unknown param: %s' % param)
						
			return es

def read_param(ls):
	assert len(ls) > 1, 'param %s requires value' % ls[0]
	return ls[1]

def read_2params(ls):
	assert len(ls) > 2, 'param %s requires two values' % ls[0]
	return ls[1], ls[2]
	
def read_3params(ls):
	assert len(ls) > 3, 'param %s requires three values' % ls[0]
	return ls[1], ls[2], ls[3]
	
	
	
	
	
