import re, os, sys, math, logging
from math3d import Vector
from ..util.units import Units
from ..structures.cell import Cell
from ..structures.atom_vector import AtomVector, AtomKeys

DFT_PARAMS = 'DFT PARAMETERS'

dft_regex = '\s+[0-9]+\s+([0-9]+)\s+([A-Z][a-z]?)\s+([0-9]{1,3}\.[0-9]{4})\s+[0-9\.\-]+'

class CrystalOut:
	OVERLAP = 0
	DENSITY = 1

	def __init__(self):
		self.valence_map = {}

			
	@staticmethod
	def from_string(datastring):
		lines = datastring.splitlines()
		co = CrystalOut()
		vm = {}
		for n in xrange(len(lines)):
			if DFT_PARAMS in lines[n]:
				break
		n += 3
		for k in xrange(n, len(lines)):
			m = re.match(dft_regex, lines[k])
			if m:
				vm[m.group(2)] = int(float(m.group(3)))
		co.valence_map = vm
		logging.debug(u'ATOM VALENCE:')
		for k in co.valence_map.keys():
			logging.debug(u'%s: %d' % (k, co.valence_map[k]))
		return co
			

	@staticmethod
	def from_file(name):
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Reading output from %s' % name)
		logging.info(u'*********************************************')
		logging.info(u'')
		with open(name) as f:
			return CrystalOut.from_string(f.read())


