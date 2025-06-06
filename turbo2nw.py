#!/usr/bin/python
import re, sys, math, logging, random
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.turbomole_format.turbo_basis import TurboBasis
from qcldm.gauss_functions.gauss_function import GaussFunction, GaussFunctionNormed, GaussFunctionContracted
from qcldm.util.log_colorizer import init_log
from qcldm.util.units import Units
from qcldm.util.elements import ELEMENTS
from qcldm.structures.atom_vector import AtomKeys
from qcldm.atom.shells import Shells

TEMPLATE = '''start cluster
title "cluster"
geometry "full-cluster" units au noautoz
symmetry c1
%%COORD%%
end


basis spherical
%%BASIS%%
end

ecp
%%ECP%%
end

so
%%SOECP%%
end

set geometry "full-cluster"
set geometry:actlist %%FRATOM%%:%%LRATOM%%
charge %%CHARGE%%
scf
 vectors input  hcore
 vectors output ./cluster.movecs
 uhf
 maxiter 100
end
task scf energy
dft
 XC pbe0
 CONVERGENCE fast
 vectors input ./cluster.movecs
 vectors output ./cluster.movecs
 iterations 200
grid lebedev 350 17 becke;tolerances accCoul 20 tol_rho 20
end

task dft gradient
'''

def dummy_basis():
	gc = GaussFunctionContracted()
	gc.fs.append([1, GaussFunctionNormed(100, 0)])
	return TurboBasis('', [gc])


def dummy_ecp(name):
	return TurboBasis.TurboEcp('', ELEMENTS[name].number, TurboBasis.EcpPart([[0, GaussFunction(1, 2)]]), [TurboBasis.EcpPart([[0.00000001, GaussFunction(1, 2)]])], [])

def build_basis(c, specs, placeholder):
	basis_str = ''
	ecp_str = ''
	so_str = ''
#	print specs
	elbasmap = {}
	elecpmap = {}
	for spec in specs:
#		print(spec)
		el = None
		bas = None
		if spec[0] == 'none' and spec[1] == None:
#			print(spec)
#			continue
			el = placeholder
			bas = dummy_basis()
			elecpmap[el] = dummy_ecp(placeholder)
		elif spec[0] == 'none':
			el = spec[1].split()[0]
			bas = dummy_basis()
		else:
			el = spec[0].split()[0]
			bas = c.bases[spec[0]]
		el = el[0].upper() + el[1:].lower()
		assert el not in list(elbasmap.keys()) or bas.name == elbasmap[el].name, 'multiple bases for %s: %s AND %s' % (el, bas.name, elbasmap[el].name)
		elbasmap[el] = bas
		if spec[1] != None:
			assert el not in list(elecpmap.keys()), 'multiple ecps for %s' % el
			elecpmap[el] = c.ecps[spec[1]]
	
	for el in sorted(elbasmap.keys()):
		bas = elbasmap[el]
		for f in bas.functions:
			l = Shells.SHELLS[f.fs[0][1].l]
			basis_str += el + ' ' * (6 - len(el)) + l + '\n'
			for c, g in f.fs:
				basis_str += "%18.10f %15.8f" % (g.a, c) + '\n'

	for el in sorted(elecpmap.keys()):
		ecp = elecpmap[el]
		ecp_str += el + ' ' * (2 - len(el)) + ' nelec ' + str(ecp.ncore) + '\n'
		ecp_str += el + ' ' * (2 - len(el)) + ' ul\n'
		for c, g in ecp.local.functions:
			ecp_str += '    %2d %15.8f %15.8f\n' % (g.l, g.a, c)
		for l, part in enumerate(ecp.semilocal):
			ecp_str += el + ' ' * (2 - len(el)) + ' ' + Shells.SHELLS[l] + '\n'
			for c, g in part.functions:
				ecp_str += '    %2d %15.8f %15.8f\n' % (g.l, g.a, c)

	for el in sorted(elecpmap.keys()):
		ecp = elecpmap[el]
		if not ecp.spinorbit:
			continue
		for l, part in enumerate(ecp.spinorbit):
			l += 1
			so_str += el + ' ' * (2 - len(el)) + ' ' + Shells.SHELLS[l] + '\n'
			for c, g in part.functions:
				so_str += '    %2d %15.8f %15.8f\n' % (g.l, g.a, c)

	return basis_str, ecp_str, so_str



init_log(sys.argv)

c = ControlFormat.from_path('.')
#for a in c.cell.atoms:
#	print a

placeholder = 'He'


with open('test.nw', 'w') as f:
	coord_block = ''
	specs = set()
	charge = 0.
	embcount = {}
	bqbases = []
	smap = c.species_map()
	first_real_atom = None
	last_real_atom = None
	for i, a in enumerate(c.cell.atoms):
		embedded = AtomKeys.ESTIMATED_CHARGE in list(a.data().keys())
		if not embedded:
			if first_real_atom is None:
				first_real_atom = str(i + 1)
			last_real_atom = str(i + 1)
		name = a.name()
		if name.lower() == 'q':
#			name = 'bq'
			name = placeholder
			a.data()[AtomKeys.ESTIMATED_CHARGE] += ELEMENTS[name].number
		if embedded:
			if name not in list(embcount.keys()):
				embcount[name] = 0
			embcount[name] = embcount[name] + 1
			if name != 'bq':
				name += str(embcount[name])
			real_basisname = smap[i + 1][0]
			if name == 'bq' and real_basisname != 'none':
				bas_to_copy = c.bases[real_basisname]
				if real_basisname not in bqbases:
					bqbases.append(real_basisname)
					bqbasisname = 'bq%d %s' % (bqbases.index(real_basisname) + 1, real_basisname.split()[1])
					bqbasis = TurboBasis(bqbasisname, bas_to_copy.functions)
					c.bases[bqbasisname] = bqbasis
					smap[i + 1] = bqbasisname, smap[i + 1][1]
				name = 'bq%d' % (bqbases.index(real_basisname) + 1)

		coord_block += "%3s %15.10f %15.10f %15.10f" % (name, a.position()[0] / Units.BOHR, a.position()[1] / Units.BOHR, a.position()[2] / Units.BOHR)

		if embedded:
			ncore = c.ecps[smap[i + 1][1]].ncore if smap[i + 1][1] else 0
			emcharge = a.data()[AtomKeys.ESTIMATED_CHARGE] + ncore
			real_basisname = smap[i + 1][0]
			if 'bq' in name and real_basisname != 'none':
				if emcharge == 0:
					emcharge = 1e-8
			charge += emcharge
			coord_block += ' charge %13.8f' % emcharge
		coord_block += "\n"
		specs.add(smap[i + 1])
	charge =  charge - round(charge)
	print(charge)	
	chargestr = "%13.8f" % charge
	basis_block, ecp_block, so_block = build_basis(c, specs, placeholder)
	f.write(TEMPLATE.replace('%%COORD%%', coord_block[:-1]).replace('%%BASIS%%', basis_block[:-1]).replace('%%ECP%%', ecp_block[:-1]).replace('%%SOECP%%', so_block[:-1]).replace('%%CHARGE%%', chargestr).replace('%%FRATOM%%', first_real_atom).replace('%%LRATOM%%', last_real_atom))























	











