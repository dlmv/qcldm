import os, logging
from ..structures.cluster_embedding import Cluster
from ..util.fileutils import make_dir
from ..structures.atom_vector import AtomVector, AtomKeys
from ..structures.cell import Cell
from ..atom.basis import Basis
from ..atom.pseudo_potential import SeparablePseudoPotential
from ..openmx_format.pao_format import PAO
from ..openmx_format.vps_format import VPS
from functools import reduce



def rewrite_float(f):
	s = 'abcdefghij'
	fs = "%5.3f" % f
	res = ''
	for c in fs:
		if c.isdigit():
			res += s[int(c)]
		elif c == '.':
			res += 'o'
		elif c == ' ':
			res += 's'
		else:
			res += 'x'
	return res

def rewrite_floats(fs):
	return reduce(lambda res, f: res + rewrite_float(f), fs, '')

def modified_basis(basis, a):
	orig_valence = a.origin.data()[AtomKeys.ESTIMATED_VALENCE]
	new_basis = Basis()

	if a.valence == orig_valence and a.charge == orig_valence:
		return basis
	elif a.valence != 0:
		new_basis.z = basis.z + a.charge - orig_valence
		new_basis.etotal = new_basis.z
		new_basis.eval = basis.eval + a.charge - orig_valence
#		new_basis.vden = basis.vden - basis.estimate_upper_vden()  * ((basis.eval - (new_basis.eval + new_basis.z - new_basis.etotal)) / basis.estimate_valence())
		new_basis.vden = basis.vden - basis.estimate_upper_vden()  * ((basis.estimate_valence() - a.valence) / basis.estimate_valence())
		new_basis.cutoff = basis.cutoff
		new_basis.functions = basis.functions
		new_basis.grid = basis.grid
		return new_basis
	else:
		new_basis.z = a.charge
		new_basis.etotal = 0
		new_basis.eval = 0
		new_basis.vden = basis.vden * 0
		new_basis.cutoff = basis.cutoff
		new_basis.functions = basis.functions
		new_basis.grid = basis.grid
		return new_basis

def modified_pp(pp, a):
	orig_valence = a.origin.data()[AtomKeys.ESTIMATED_VALENCE]
	new_pp = SeparablePseudoPotential()

	if a.valence == orig_valence and a.charge == orig_valence:
		return pp
	elif a.valence != 0:
		new_pp.z = pp.z + a.charge - orig_valence
		new_pp.etotal = new_pp.z
		new_pp.eval = pp.eval + a.charge - orig_valence
		new_pp.pcc = pp.pcc
		new_pp.grid = pp.grid
		new_pp.loc = pp.loc * (new_pp.eval / pp.eval)
		new_pp.blochl_projectors = pp.blochl_projectors
		return new_pp
	else:
		new_pp.z = a.charge
		new_pp.etotal = 0
		new_pp.eval = 0
		new_pp.pcc = None
		new_pp.grid = pp.grid
		new_pp.loc = pp.loc * (a.charge / pp.eval)
		new_pp.blochl_projectors = {}
		return new_pp

def rewrite_dat_file(dat, cell, newspecies, datadir):
	dat.data_path = "."
	dat.species = newspecies
	dat.cell = cell
	dat.system_name = 'cluster'


	dat.to_file(os.path.join(datadir, "cluster.dat"))


def rewrite_files(dat, cluster, datadir):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Creating cluster input')
	logging.info('*********************************************')
	logging.info('')
	newspecies = {}
	newatoms = []
	paodir = os.path.join(datadir, "PAO")
	vpsdir = os.path.join(datadir, "VPS")
	make_dir(paodir)
	make_dir(vpsdir)
	for i, a in enumerate(cluster.atoms):
		sp = dat.species[a.origin.name()]
#		fstring = rewrite_float(a.charge)
		fstring = rewrite_floats([a.valence, a.charge])
		new_basis = modified_basis(sp.basis, a)
		PAO.to_file(new_basis, os.path.join(paodir, '%s.pao' % (fstring + sp.pao)))
		new_pp = modified_pp(sp.pp, a)
		VPS.to_file(new_pp, os.path.join(vpsdir, '%s.vps' % (fstring + sp.vps)))
		newsp = sp.prefixed(fstring, datadir)
		newspecies[newsp.name] = newsp
		a = AtomVector(newsp.name, a.origin.position())
		newatoms.append(a)
		logging.debug('  %d/%d' % (i + 1, len(cluster.atoms)))
	rewrite_dat_file(dat, Cell(newatoms, []), newspecies, datadir)
