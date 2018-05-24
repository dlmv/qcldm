import re, os, sys, logging
from ..structures.atom_vector import AtomKeys

def mulliken_population(a, dm, olp):
	a = a.relative()
	numorbs = a.data()[AtomKeys.ORBITAL_COUNT]
	pop = 0
	keys = set(dm['a']['a']['re'].keys()) & set(olp.keys())
	for n in range(numorbs):
		orbpop = 0
		orb = (a.tuple_data(), n + 1)
		for key in keys:
			if orb in key:
				for tdm in dm['a']['a']['re'], dm['b']['b']['re']:
					orbpop += tdm[key] * olp[key]
		pop += orbpop
	return pop

def mulliken_overlap(a1, a2, dm, olp):
	a2 = a2.relative(a1)
	a1 = a1.relative()
	numorbs1 = a1.data()[AtomKeys.ORBITAL_COUNT]
	numorbs2 = a2.data()[AtomKeys.ORBITAL_COUNT]
	pop = 0
	keys = set(dm['a']['a']['re'].keys()) & set(olp.keys())
	for n1 in range(numorbs1):
		orb1 = (a1.tuple_data(), n1 + 1)
		for n2 in range(numorbs2):
			orb2 = (a2.tuple_data(), n2 + 1)
			key = tuple(sorted([orb1, orb2]))
			if key in keys:
				for tdm in dm['a']['a']['re'], dm['b']['b']['re']:
					pop += tdm[key] * olp[key]
#	print a1, a2, pop
	return pop


#def mayer_index(a1, a2, dms, olp):
#	a2 = a2.relative(a1)
#	a1 = a1.relative()
#	numorbs1 = a1.data()[AtomKeys.ORBITAL_COUNT]
#	numorbs2 = a2.data()[AtomKeys.ORBITAL_COUNT]
#	keys = set(dms[0].keys()) & set(olp.keys())
#	keys1 = set()
#	keys2 = set()
#	a2a = a2.relative()
#	keys2a = set()
#	for k in keys:
#		if (a1.tuple_data() == k[0][0] or a1.tuple_data() == k[1][0]):
#			keys1.add(k)
#		if (a2.tuple_data() == k[0][0] or a2.tuple_data() == k[1][0]):
#			keys2.add(k)
#		if (a2a.tuple_data() == k[0][0] or a2a.tuple_data() == k[1][0]):
#			keys2a.add(k)
#	pop = 0
#	for n1 in range(numorbs1):
#		orb1 = (a1.tuple_data(), n1 + 1)
#		for n2 in range(numorbs2):
#			ps12 = 0
#			orb2 = (a2.tuple_data(), n2 + 1)
#			orb2a = (a2a.tuple_data(), n2 + 1)
#			for key1 in keys1:
#				if orb1 in key1:
#					orbk = key1[1] if orb1 == key1[0] else key1[0]
#					key2 = tuple(sorted([orbk, orb2]))
#					ak = orbk[0]
#					aka = (ak[0], ak[1]-a2.na, ak[2]-a2.nb, ak[3]-a2.nc)
#					orbka = (aka, orbk[1])
#					key2a = tuple(sorted([orbka, orb2a]))
#					for dm in dms:
#						if key2 in keys2:
#							ps12 += dm[key1] * olp[key2] + dm[key2] * olp[key1]
#						elif a2 != a2a and key2a in keys2a:
#							ps12 += dm[key1] * olp[key2a] + dm[key2a] * olp[key1]
#			pop += ps12**2
#	pop /= 2.
#	return pop



