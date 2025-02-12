import sys, os, re
import numpy as np
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.core.periodic_table import DummySpecies

def read_crystal_cell(filename):
	with open(filename) as f:
		natoms = -1
		xyz = []
		for line in f:
			if 'LATTICE PARAMETERS (ANGSTROMS AND DEGREES)' in line:
				break
		for i in range(3):
			line = next(f)
		param_list = [float(x) for x in line.strip().split()]
		vectors = param_list[0:3]
		angles = param_list[3:6]
		for line in f:
			if 'ATOMS IN THE UNIT CELL' in line:
				natoms = int(line.strip().split()[-1])
				break
		for i in range(2):
			next(f)
		for i in range(natoms):
			ls = next(f).strip().split()
			atom_name = 'a' + ls[2]
			atom_name = atom_name[0:1].upper() + atom_name[1:].lower()
			coords = [float(x) for x in ls[4:7]]
			xyz.append([atom_name] + coords)
	return vectors, angles, xyz

		
def prepare_input(atoms, group_class, group, vectors, angles, d12_template, d12_out):
	text = str(group) + '\n'
	print(group_class, group)
	if group_class == 'monoclinic':
		text += '    %f %f %f %f\n' % (vectors[0], vectors[1], vectors[2], angles[1])
	elif group_class == 'triclinic':
		text += '    %f %f %f %f %f %f\n' % (vectors[0], vectors[1], vectors[2], angles[0], angles[1], angles[2])
	elif group_class == 'orthorhombic':
		text += '    %f %f %f\n' % (vectors[0], vectors[1], vectors[2])
	else:
		assert False, group_class
	text += str(len(atoms)) + '\n'
	for name, x, y, z in atoms:
		text += ' %3s %14.8f %14.8f %14.8f\n' % (name, x,y,z)
	new_input = None
	with open(d12_template) as fin:
		template = fin.read()
		m = re.search('\d+\s*\n\s*([\d.]+\s*){3,6}\s*\n\s*\d+\n(\s*\d+(\s+[\d.\-+E]+){3})+', template)
		new_input = re.sub('\d+\s*\n\s*([\d.]+\s*){3,6}\s*\n\s*\d+\n(\s*\d+(\s+[\d.\-+E]+){3})+\n', text, template, 1)
		assert new_input != template
	with open(d12_out, 'w') as fout:
		fout.write(new_input)

crystal_out = sys.argv[1]
crystal_d12 = sys.argv[1].replace('.out', '.d12')
assert crystal_out != crystal_d12

new_d12 = sys.argv[1].replace('.out', '_sym.d12')

vectors, angles, xyz = read_crystal_cell(crystal_out)
a,b,c = vectors
al,be,ga = angles
lattice = Lattice.from_parameters(a,b,c,al,be,ga)
species = [x[0] for x in xyz]
coords = [x[1:] for x in xyz]
structure = Structure(lattice, species, coords)
sga = SpacegroupAnalyzer(structure)
sga = SpacegroupAnalyzer(sga.get_refined_structure())
group = sga.get_space_group_number()
group_class = sga.get_crystal_system()
ss = sga.get_symmetrized_structure()
#print(ss)
vectors = [ss.lattice.a, ss.lattice.b, ss.lattice.c]
angles = [ss.lattice.alpha, ss.lattice.beta, ss.lattice.gamma]
atoms = []
symmpos = sga.get_symmetry_operations()
#for so in symmpos:
#	print('============================')
#	print(so.rotation_matrix)
#	print(so.translation_vector)
for i in ss.equivalent_indices:
	site = ss.sites[i[0]]
	x,y,z = site.frac_coords
#	print(site.as_dict())
	name = str(int(site.specie.oxi_state))
#	print(name)
#	name, x,y,z = xyz[i[0]]
	atoms.append([name, x,y,z])


prepare_input(atoms, group_class, group, vectors, angles, crystal_d12, new_d12)


















































