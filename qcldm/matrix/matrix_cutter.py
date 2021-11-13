import os, logging
from ..util.units import Units
from ..structures.atom_vector import AtomKeys
from ..util.fileutils import make_dir
from ..util.xyz_format import write_xyz

def write_reduced(dm, num, atoms, rc, dirname):
	make_dir(dirname)
	cluster = cut_cluster(num, atoms, rc)
	write_dm(dm, cluster, os.path.join(dirname, "DM.out"))
	write_xyz(atoms, os.path.join(dirname, "full.xyz"))
	write_xyz(cluster, os.path.join(dirname, "cluster.xyz"))

def cut_cluster(num, atoms, rc):
	rc /= Units.UNIT / Units.BOHR
	ca = atoms[0].cell.cell[num - 1]
	res = [ca]
	for a in atoms:
		if a == ca:
			continue
		r = a.distance(ca)
		ra = a.data()[AtomKeys.CUTOFF]
		maxr = ra + rc
		if r < maxr:
			res.append(a)
	return res

def write_dm(dm, cluster, filename):
	with open(filename, 'w') as f:
		f.write(str(len(cluster)) + "\n")
		for a in cluster:
			k = Units.UNIT / Units.BOHR
			f.write("%3d %3s %12.8f %12.8f %12.8f\n" % (a.num, a.name(), a.position().x * k, a.position().y * k, a.position().z * k))
		for i1 in ["a", "b"]:
			for i2 in ["a", "b"]:
				if i2 not in dm[i1].keys():
					continue
				for c in ["re", "im"]:
					if c not in dm[i1][i2].keys():
						continue
					f.write("\n%s Density matrix spin = %s%s" % (c, i1, i2))
					tdm = dm[i1][i2][c] if i1 != 'b' else dm[i2][i1][c]
					negative = (i1 == 'b' and i2 == 'a' and c == 'im')
					write_matrix(tdm, cluster, f, negative)

import traceback

def write_matrix(dm, cluster, f, negative=False):
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Writing reduced matrix')
	logging.info(u'*********************************************')
	logging.info(u'')
	import traceback

	header = '\n************'
	for i, a in enumerate(cluster):
		n = a.data()[AtomKeys.ORBITAL_COUNT]
		for o in range(n):
			header += " %13d %3d" % (i + 1, o + 1)
	f.write(header)
	f.write("\n")
	sp = 0
	for i, a1 in enumerate(cluster):
		a10 = a1.relative()
		n1 = a1.data()[AtomKeys.ORBITAL_COUNT]
		for o1 in range(n1):
			f.write(" %13d %3d" % (i + 1, o1 + 1))
			for j, a2 in enumerate(cluster):
				a20 = a2.relative(a1)
				n2 = a2.data()[AtomKeys.ORBITAL_COUNT]
				for o2 in range(n2):
					key = tuple(sorted([(a10.tuple_data(), o1 + 1), (a20.tuple_data(), o2 + 1)]))
					el = dm.get(key) or 0
					if negative:
						el *= -1
					f.write(" %17.12f" % el)
					if key[0] == key[1]:
						sp += el
			f.write("\n")
		logging.debug(u'Progress: %d/%d' % (i + 1, len(cluster)))








