import os
class EmbeddingSettings:
	def __init__(self):
		self.center = -1
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
						es.center = int(read_param(ls))
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
					elif param == 'basisfull':
						atom, basis = read_2params(ls)
						es.basis_map[(atom, True)] = basis
					elif param == 'basisborder':
						atom, basis = read_2params(ls)
						es.basis_map[(atom, False)] = basis
					elif param == 'ecpfull':
						atom, ecp = read_2params(ls)
						es.ecp_map[(atom, True)] = ecp
					elif param == 'ecpborder':
						atom, ecp = read_2params(ls)
						es.ecp_map[(atom, False)] = ecp
						
						
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
	
	
	
	
	
