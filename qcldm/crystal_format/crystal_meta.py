import re, os, sys, math, logging

from ..util.units import Units

dm_regex = 'BASISSET\s+1\s+1\s+59\s+256\s+END'
olp_regex = 'BASISSET\s+1\s+1\s+60\s+256\s+END'


class CrystalMeta():
	def __init__(self):
		self.out_file = None
		self.dm_file = None
		self.olp_file = None

	def load(self, dr):
		logging.info('')
		logging.info('*********************************************')
		logging.info('  Looking for crystal output files')
		logging.info('*********************************************')
		logging.info('')
		files = os.listdir('.')
		outs = []
		for o in sorted([x for x in files if x.endswith('.out')]):
			if os.path.exists(o[:-4] + '.d12'):
				outs.append(o)
		if len(outs) == 1:
			self.out_file = outs[0]
			logging.info('  Out found: %s' % self.out_file)
		elif len(outs) == 0:
			logging.error('  No outs found!!')
			assert False
		elif len(outs) > 1:
			logging.warn('  Multiple outs found!!')
			logging.error('  Remove some, or fix my code!')
			assert False
		outps = sorted([x for x in files if x.endswith(self.out_file + 'p')])
		for outp in outps:
			d3inp = outp[:-len(self.out_file) - 2] + '.d3'
			with open(d3inp) as f:
				data = f.read()
#				head = ''
#				lim = 20
#				for line in f.xreadlines():
#					lim -=1
#					if lim < 0:
#						break
#					head += line
				if re.search(dm_regex, data):
					self.dm_file = outp
					logging.info('  DM found: %s' % self.dm_file)

				elif re.search(olp_regex, data):
					self.olp_file = outp
					logging.info('  OLP found: %s' % self.olp_file)

					
				
			
		


