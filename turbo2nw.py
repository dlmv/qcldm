#!/usr/bin/python
import re, sys, math, logging, random
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.turbomole_format.turbo_basis import TurboBasis
from qcldm.gauss_functions.gauss_function import GaussFunctionNormed, GaussFunctionContracted
from qcldm.util.log_colorizer import init_log
from qcldm.util.units import Units
from qcldm.structures.atom_vector import AtomKeys
from qcldm.atom.shells import Shells

TEMPLATE = '''start cluster 
title "cluster" 
geometry "full-cluster" units au 
%%COORD%%
end  


basis spherical
%%BASIS%%
end
 
ecp 
%%ECP%%
end

set geometry "full-cluster"
charge %%CHARGE%%
dft  
 XC pbe0
 vectors input hcore
 vectors output ./cluster.movecs
 iterations 200
grid xfine
tolerances accCoul 20 tol_rho 20
end  

task dft
'''

def dummy_basis():
	gc = GaussFunctionContracted()
	gc.fs.append([1, GaussFunctionNormed(100000, 0)])
	return TurboBasis('', [gc])
	

def build_basis(c, specs):
	basis_str = ''
	ecp_str = ''
#	print specs
	elbasmap = {}
	elecpmap = {}
	for spec in specs:
#		print spec
		if spec[0] == 'none' and spec[1] == None:
			continue
		el = None
		bas = None
		if spec[0] == 'none':
			el = spec[1].split()[0]
			bas = dummy_basis()
		else:
			el = spec[0].split()[0]
			bas = c.bases[spec[0]]
		el = el[0].upper() + el[1:].lower()
		assert el not in elbasmap.keys(), 'multiple bases for %s: %s AND ' % (el, bas.name, elbasmap[el])
		elbasmap[el] = bas
		if spec[1] != None:
			assert el not in elecpmap.keys(), 'multiple ecps for %s' % el
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

	return basis_str, ecp_str
		
#	for elm in sorted(elms):
#		elm = elm[0].upper() + elm[1:].lower()
#		found = False
#		for bn in c.bases.keys():
#			atomname = bn.split()[0]
#			atomname = atomname[0].upper() + atomname[1:].lower() 
#			if atomname == elm:
#				assert found, 'multiple basis set for '
#		b = c.bases[bn]
#		for f in b.functions:
#			print f




init_log(sys.argv)

c = ControlFormat.from_file('control')
#for a in c.cell.atoms:
#	print a



with open('test.nw', 'w') as f:
	coord_block = ''
	specs = set()
	charge = 0.
	embcount = {}
	for i, a in enumerate(c.cell.atoms):
		embedded = AtomKeys.ESTIMATED_CHARGE in a.data().keys()
		name = a.name()
		if name.lower() == 'q':
			name = 'bq'
		if embedded:
			if name not in embcount.keys():
				embcount[name] = 0
			embcount[name] = embcount[name] + 1
			if name != 'bq':
				name += str(embcount[name])
		coord_block += "%3s %15.10f %15.10f %15.10f" % (name, a.position()[0] / Units.BOHR, a.position()[1] / Units.BOHR, a.position()[2] / Units.BOHR)
		if embedded:
			ncore = c.ecps[c.species_map()[i + 1][1]].ncore if c.species_map()[i + 1][1] else 0
			emcharge = a.data()[AtomKeys.ESTIMATED_CHARGE] + ncore
			charge += emcharge
			coord_block += ' charge %13.8f' % emcharge
		coord_block += "\n"
		specs.add(c.species_map()[i + 1])
	charge =  charge - round(charge)
	print charge	
	chargestr = "%13.8f" % charge
	basis_block, ecp_block = build_basis(c, specs)
	f.write(TEMPLATE.replace('%%COORD%%', coord_block[:-1]).replace('%%BASIS%%', basis_block[:-1]).replace('%%ECP%%', ecp_block[:-1]).replace('%%CHARGE%%', chargestr))























	











