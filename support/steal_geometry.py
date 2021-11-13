#!/usr/bin/python 
import os, sys, re, numpy as np



path1 = '.'
path2 = sys.argv[1]

class Atom:
	def __init__(self, name, coords):
		self.coords = coords
		self.name = name
		self.emb_name = None
		self.emb_charge = 0
		self.emb_min = 0
		self.emb_max = 0
		self.emb_group = None
		self.basis = None
		self.ecp = None
	
	def real_name(self):
		return self.emb_name if self.emb_name else self.name

	def copy(self):
		a = Atom(self.name, self.coords)
		a.emb_name = self.emb_name
		a.emb_charge = self.emb_charge
		a.emb_min = self.emb_min
		a.emb_max = self.emb_max
		a.emb_group = self.emb_group
		a.basis = self.basis
		a.ecp = self.ecp
		return a
		
def distance(a1, a2):
	assert len(a1.coords) == 3
	assert len(a2.coords) == 3
	return sum([(a - b) ** 2 for a, b in zip(a1.coords, a2.coords)])**0.5

def steal_coords(as1, as2):
	deleted_atoms = []
	for i1, a1 in enumerate(as1):
		row = [0] * len(as1)
		i2, a2 = min(enumerate(as2), key=lambda x: distance(a1, x[1]))
		if a1.name != a2.name or a1.emb_name != a2.emb_name:
			if a1.emb_name != a2.name:
				if a1.emb_name != 'q':
					if a2.emb_name == 'q':
						copy = a1.copy()
						deleted_atoms.append(copy)
						print 'deleted %d %s %s' % (i1, a1.name, a1.emb_name)
						a1.name = 'zz'
						a1.emb_name = 'q'
						a1.ecp = None
						a1.emb_charge = 0
						a1.emb_min = 0
						a1.emb_max = 0
					else:
						assert False

		a1.coords = a2.coords
	return deleted_atoms


def read_coord_and_embedding(path):
	atoms = []
	with open(os.path.join(path, 'coord')) as coord:
		lines = coord.read().splitlines()
		assert lines[0] == '$coord'
		assert lines[-1] == '$end'
		for line in lines[1:-1]:
			ls = line.split()
			x, y, z = [float(r) for r in ls[:3]]
			name = ls[3]
			a = Atom(name.lower(), [x,y,z])
			atoms.append(a)
		
	with open(os.path.join(path, 'embedding')) as embedding:
		lines = embedding.read().splitlines()
		pos = 0
		for line in lines:
			ls = line.split()
			name = ls[0]
			charge = float(ls[1])
			while atoms[pos].name.lower() != 'zz':
				pos += 1
			
			atoms[pos].emb_name = name.lower()
			atoms[pos].emb_charge = charge
			pos += 1
	if os.path.exists(os.path.join(path, 'embedding.restart')):
		with open(os.path.join(path, 'embedding.restart')) as embedding:
			lines = embedding.read().splitlines()
			pos = 0
			for line in lines:
				ls = line.split()
				name = ls[0]
				charge = float(ls[1])
				
				while atoms[pos].name.lower() != 'zz':
					pos += 1
				
				atoms[pos].emb_min = float(ls[2])
				atoms[pos].emb_max = float(ls[3])
				atoms[pos].emb_group = ls[4]
				pos += 1
	return atoms


def read_control_new(path, atoms):
	with open(os.path.join(path, 'control')) as control:
		lines = control.read().splitlines()
		reading_atoms = False
		name = None
		cur_list = []
		new_lines = []
		for line in lines:
			if '$atoms' in line:
				reading_atoms = True
				new_lines.append(line)
				continue
			elif '$' in line:
				reading_atoms = False
				new_lines.append(line)
				continue
			if reading_atoms:
				m = re.match("([a-z]+)\\s+([\d\\-, ]+?) +\\\\", line)
				if m:
					name = m.group(1)
					ls = m.group(2).strip().split(',')
					cur_list = []
					for part in ls:
						m = re.match('(\d+)-(\d+)', part)
						if m:
							cur_list.extend(range(int(m.group(1)), int(m.group(2)) + 1))
							continue
						m = re.match('(\d+)', part)
						if m:
							cur_list.append(int(m.group(1)))
							continue
						assert True, 'AAAAAAAAAAAAAAA'
					for i in cur_list:
						assert name == atoms[i-1].real_name()
				else:
					m = re.match(" *(basis|ecp) *= *([a-z]{1,2} ?[^ ]+)( +\\\\)?", line)
					if m:
						if m.group(1) == 'basis':
							for i in cur_list:
								atoms[i-1].basis = m.group(2)
						elif m.group(1) == 'ecp':
							for i in cur_list:
								atoms[i-1].ecp = m.group(2)
			else:
				new_lines.append(line)
				continue
		return new_lines
			
				
			
def numstr(numbers):
	lastfirst = numbers[0]
	current = numbers[0]
	res = ''
	for x in numbers[1:] + [999]:
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


def write_control_new(path, atoms, lines):
	try:
		os.makedirs(path)
	except Exception:
		pass
	species_map = {}
	for i, a in enumerate(atoms):
		specie = "#".join([a.real_name(), a.basis, a.ecp] if a.ecp else [a.real_name(), a.basis])
		if specie not in species_map:
			species_map[specie] = []
		species_map[specie].append(i + 1)
	atompart = ''
	for specie in sorted(species_map.keys(), key=lambda x: species_map[x][0]):
		ls = specie.split('#')
		name, basis, ecp = ls if len(ls) == 3 else ls + [None]
		line = '{:<3}'.format(name) + numstr(species_map[specie])
		line = '{:<79}\\'.format(line)
		atompart += line + '\n'
		line = '   basis =' + basis
		if ecp:
			line = '{:<79}\\'.format(line)
		atompart += line + '\n'
		if ecp:
			line = '   ecp   =' + ecp
			atompart += line + '\n'
	with open(os.path.join(path, 'control'), 'w') as control:
		for line in lines:
			control.write(line + '\n')
			if '$atoms' in line:
				control.write(atompart)
		
		
			

def write_coord_and_embedding(path, atoms, atomsd):
	try:
		os.makedirs(path)
	except Exception:
		pass
	with open(os.path.join(path, 'coord'), 'w') as coord:
		coord.write('$coord\n')
		for a in atoms:
			coord.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(a.coords[0], a.coords[1], a.coords[2], a.name))
		coord.write('$end\n')
	with open(os.path.join(path, 'coord_ca'), 'w') as coord:
		coord.write('$coord\n')
		for a in atoms:
			coord.write("  {:15.10f}  {:15.10f}  {:15.10f}  {:3}\n".format(a.coords[0], a.coords[1], a.coords[2], a.real_name()))
		coord.write('$end\n')
	with open(os.path.join(path, 'embedding'), 'w') as embedding:
		for a in atoms + atomsd:
			if a.emb_name:
				embedding.write("{:3}  {:9.5f}\n".format(a.emb_name, a.emb_charge))
	with open(os.path.join(path, 'embedding.start'), 'w') as embedding:
		for a in atoms + atomsd:
			if a.emb_name:
				embedding.write("{:3}  {:9.5f}  {:9.5f}  {:9.5f} {}\n".format(a.emb_name, a.emb_charge, a.emb_min, a.emb_max, a.emb_group))
				
				
atoms = read_coord_and_embedding(path1)
atoms2 = read_coord_and_embedding(path2)

newlines = read_control_new(path1, atoms)


atomsd = steal_coords(atoms, atoms2)

n = 1
for a in atoms:
	if a.emb_group:
		a.emb_group = "SINGLE_%d" % n
		n += 1
for a in atomsd:
	if a.emb_group:
		a.emb_group = "SINGLE_%d" % n
		n += 1

dirname = os.path.basename(os.path.dirname(path2))
write_coord_and_embedding('stolen_from_%s' % dirname, atoms, atomsd)
write_control_new('stolen_from_%s' % dirname, atoms, newlines)





















