
import re, os, logging, shutil
import configparser

from ..matrix.matrix_cutter import write_reduced
from ..gauss_functions.gauss_formats import GaussFormat

OPTIONS = '''
ObjVersion = 1
Libdir = %LIB%
OutputName = .\\%OUT%
LogDir = .\\log\\
NumberOfAlphaElectrons = 0
NumberOfBetaElectrons = 0
NumberOfUAlphaElectrons = 0
NumberOfUBetaElectrons = 0
CalcCoreProps = -1
SaveCoreInts = 0
CalculateDipole = 0
CalculateMoments = 0
SaveValenceInts = 0
AnalyzeOrbitals = -1
AtomCoreNum = 0
SecondAtomNumber = 1
UseDMcompl = 0
Rc = %RC%
RcStep = 0.1
RcStepNumber = 1
Is4c = 0
CalcNonDiag = -1
CalcDiag4c = 0
IncludeHFDens = -1
nFroz = 0
nUFroz = 0
nFrozVirt = 0
nUFrozVirt = 0
Smin = 1E-7
Sq = 1.4
Smax = 1000000
Pmin = 1E-7
Pq = 1.4
Pmax = 1000000
Dmin = 1E-7
Dq = 1.4
Dmax = 1000000
Fmin = 1E-7
Fq = 1.4
Fmax = 1000000
RcAlgoritm = 0
RcAccuracy = 300
PrintExternal = 0
'''[1:]

def get_dirs():
	gaussdir = '.'
	libdir = '.\\lib\\'

	config = configparser.ConfigParser()
	config.read(os.path.join(os.path.expanduser("~"), ".oneprop.ini"))
	try:
		gaussdir = config.get("dirs", "gaussdir")
	except Exception:
		pass
	try:
		libdir = config.get("dirs", "libdir")
	except Exception:
		pass

	return gaussdir, libdir

def prepare_oneprop_openmx(dat, atoms, dms, num, rc):
	gaussdir, libdir = get_dirs()

	dirname = "reduced_{}_{}".format(num, rc)
	write_reduced(dms, num, atoms, rc, dirname)
	outname = '%s.out' % dat.system_name
	shutil.copy(outname, dirname)
	options = OPTIONS.replace('%RC%', str(rc)).replace('%LIB%', libdir).replace('%OUT%', outname)
	with open(os.path.join(dirname, 'Options.ini'), 'w') as f:
		f.write(options)

	basis = ''
	for s in list(dat.species.values()):
		gfile = os.path.join(gaussdir, s.pao + '.gauss')
		title, fs = read_gauss_basis(gfile)
		title = title.replace('X', s.name)
		bs = basis_to_string(fs, s.numorbs)
		basis += title
		basis += bs
		basis += '****\n'
	with open(os.path.join(dirname, 'basis.L'), 'w') as f:
		f.write(basis)

def prepare_oneprop_crystal(co, atoms, dms, num, rc):
	gaussdir, libdir = get_dirs()
	dirname = "reduced_{}_{}".format(num, rc)
	write_reduced(dms, num, atoms, rc, dirname)
	outname = co.name
	shutil.copy(outname, dirname)
	options = OPTIONS.replace('%RC%', str(rc)).replace('%LIB%', libdir).replace('%OUT%', outname)
	with open(os.path.join(dirname, 'Options.ini'), 'w') as f:
		f.write(options)
		
	basis = ''
	for k in list(co.basis.keys()):
		basis += '%s 0\n' % k
		basis += GaussFormat.basis_to_gaussian94(co.basis[k])
		basis += '****\n'
	with open(os.path.join(dirname, 'basis.L'), 'w') as f:
		f.write(basis)
		
def read_gauss_basis(gfile):
	lines = open(gfile).read().splitlines()
	title = lines[0] + '\n'
	n = 1
	cur = ''
	fs = {}
	while n < len(lines):
		line = lines[n]
		if line[0] != ' ':
			if cur:
				l = cur[0].lower()
				if l not in list(fs.keys()):
					fs[l] = []
				fs[l].append(cur)
			cur = line + '\n'
		else:
			cur += line + '\n'
		n += 1
		if n == len(lines):
			if cur:
				l = cur[0].lower()
				if l not in list(fs.keys()):
					fs[l] = []
				fs[l].append(cur)
	return title, fs

def basis_to_string(fs, nums):
	shells = "spdfgh"
	res = ''
	for l in shells:
		if l in list(fs.keys()):
			if l not in list(nums.keys()):
				continue
			for n in range(nums[l]):
				lf = fs[l][n]
				res += lf
	return res





























