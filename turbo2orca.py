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

TEMPLATE = '''!PBE0 TightSCF RIJCOSX AutoAux
%method
  IntAccX 5,5,5
  GridX 3,3,4
end
%maxcore 45000
%scf
guess hcore
maxiter 200
end
%Pal nprocs 28 end
%geom
 Constraints
  {C %%NS%% C}
 end
end
%basis
%%BASIS%%
end
*xyz %%CHARGE%% 1
%%COORDS%%
*
'''

def build_basis(c, specs):
	basis_str = ''
	ecp_str = ''
#	print specs
	elbasmap = {}
	elecpmap = {}
	for spec in specs:
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
		basis_str += ' NewGTO %d\n' % ELEMENTS[el].number
		for f in bas.functions:
			l = Shells.SHELLS[f.fs[0][1].l]
			basis_str += '      %s %d\n' % (l, len(f.fs))
			for i, (c, g) in enumerate(f.fs):
				basis_str += '      %d %18.10f %18.10f\n'  % (i+1, g.a, c)
		basis_str += ' end\n'

	for el in sorted(elecpmap.keys()):
		ecp = elecpmap[el]
		ecp_str += ' NewECP %d\n' % ELEMENTS[el].number
		ecp_str += '   N_core %d\n' % ecp.ncore
		ecp_str += '      lmax %s\n' % Shells.SHELLS[len(ecp.semilocal)].lower()
		for l, part in enumerate(ecp.semilocal + [ecp.local]):
			ecp_str += '      %s %d\n' % (Shells.SHELLS[l].lower(), len(part.functions))
			for i, (c, g) in enumerate(part.functions):
				ecp_str += '     %2d %18.10f %18.10f %d\n' % (i+1, g.a, c, g.l)
		ecp_str += ' end\n'
	return basis_str[:-1], ecp_str[:-1]



init_log(sys.argv)

el_replacements = {}
for el in 'He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn':
	el1 = ELEMENTS[ELEMENTS[el].number + 1].symbol
	el_replacements[el] = el1

for arg in sys.argv[1:]:
	m = re.match('^([A-Z][a-z]?)=([A-Z][a-z]?)$', arg)
	el_replacements[m.group(1)] = m.group(2)


c = ControlFormat.from_path('.')
#for a in c.cell.atoms:
#	print a

def numstr(numbers):
	lastfirst = numbers[0]
	current = numbers[0]
	res = ''
	for x in numbers[1:] + [9999]:
		if current == x-1:
			current = x
		elif current == lastfirst:
			res += "%d," % current
			current = x
			lastfirst = x
		else:
			res += "%d:%d," % (lastfirst, current)
			lastfirst = x
			current = x
	return res[:-1]

with open('test.inp', 'w') as f:
	coord_block = ''
	specs = set()
	main_cluster_electrons = 0
	smap = c.species_map()
	sp_replacements = {}
	emb_numbers = []
	for i, a in enumerate(c.cell.atoms):
		embedded = AtomKeys.ESTIMATED_CHARGE in list(a.data().keys())
		name = a.name()
#		print(name)
		if name.lower() == 'q':
			emb_numbers.append(i)
			coord_block += "%2s %13.8f %15.10f %15.10f %15.10f\n" % (name, a.data()[AtomKeys.ESTIMATED_CHARGE], a.position()[0] / Units.ANGSTROM, a.position()[1] / Units.ANGSTROM, a.position()[2] / Units.ANGSTROM)
		else:
			
			if embedded:
				if name in el_replacements:
					delta_val = Shells.estimate_valence_byname(el_replacements[name]) - Shells.estimate_valence_byname(name)
					name = el_replacements[name]
					if smap[i + 1] not in sp_replacements.keys():
						repl = (name.lower() + ' ' + smap[i + 1][0].split()[1], name.lower() + ' ' + smap[i + 1][1].split()[1])
						sp_replacements[smap[i + 1]] = repl
						specs.add(repl)
						c.bases[repl[0]] = c.bases[smap[i + 1][0]]
						c.ecps[repl[1]] = c.ecps[smap[i + 1][1]]
					main_cluster_electrons += delta_val
				else:
					specs.add(smap[i + 1])
				emb_numbers.append(i)
				coord_block += "%2s %15.10f %15.10f %15.10f Z = %13.8f\n" % (name, a.position()[0] / Units.ANGSTROM, a.position()[1] / Units.ANGSTROM, a.position()[2] / Units.ANGSTROM, a.data()[AtomKeys.ESTIMATED_CHARGE])
				
			else:
				coord_block += "%2s %15.10f %15.10f %15.10f\n" % (name, a.position()[0] / Units.ANGSTROM, a.position()[1] / Units.ANGSTROM, a.position()[2] / Units.ANGSTROM)
				specs.add(smap[i + 1])
				main_cluster_electrons += a.data()[AtomKeys.FULL_VALENCE]
	chargestr = str(main_cluster_electrons - c.nocc)
	basis_block, ecp_block = build_basis(c, specs)
	constr = numstr(emb_numbers)
	f.write(TEMPLATE.replace('%%COORDS%%', coord_block[:-1]).replace('%%BASIS%%', basis_block + '\n' + ecp_block).replace('%%CHARGE%%', chargestr).replace('%%NS%%', constr))























	











