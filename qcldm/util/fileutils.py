
import os, logging, errno

def ensure_file(f):
	if not os.path.isfile(f):
		logging.error(u'File not found: %s' % f)
		raise RuntimeError()
		
def ensure_files(inps, outps):
	idate = 0.
	odate = float("inf")
	for f in inps:
		if not os.path.isfile(f):
			logging.error(u'File not found: %s' % f)
			raise RuntimeError()
		idate = max(idate, os.path.getmtime(f))
	for f in outps:
		if os.path.isfile(f):
			odate = min(odate, os.path.getmtime(f))
		else:
			return True
	return idate > odate
	
def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
	
