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
		logging.info(u'')
		logging.info(u'*********************************************')
		logging.info(u'  Looking for crystal output files')
		logging.info(u'*********************************************')
		logging.info(u'')
		files = os.listdir('.')
		outs = sorted(filter(lambda x: x.endswith('.out'), files))
		if len(outs) == 1:
			self.out_file = outs[0]
			logging.info(u'  Out found: %s' % self.out_file)
		elif len(outs) == 0:
			logging.error(u'  No outs found!!')
			assert False
		elif len(outs) > 1:
			logging.warn(u'  Multiple outs found!!')
			logging.error(u'  Remove some, or fix my code!')
			assert False
		outps = sorted(filter(lambda x: x.endswith(self.out_file + 'p'), files))
		for outp in outps:
			with open(outp) as f:
				head = ''
				lim = 20
				for line in f.xreadlines():
					lim -=1
					if lim < 0:
						break
					head += line
				if re.search(dm_regex, head):
					self.dm_file = outp
					logging.info(u'  DM found: %s' % self.dm_file)

				elif re.search(olp_regex, head):
					self.olp_file = outp
					logging.info(u'  OLP found: %s' % self.olp_file)

					
				
			
		


