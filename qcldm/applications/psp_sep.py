#!/usr/bin/python
import re, sys, os, logging
sys.dont_write_bytecode = True

from ..atom.basis import Basis
from ..record_format.psp_format import read_psp
from ..record_format.func_formats import read_hfj_dft
from .psp_preparer import BASIS_DIR, PROJ_DIR

PROJ_E_LIMIT = 0.3

def get_hfj_depth(inp):
	with open(inp) as f:
		lines = f.readlines()
		for l in lines:
			if "RWL" in l:
				return float(re.split("\s+", l.split("=")[1].strip())[1])
	return 0

def get_hfj_energies(inp):
	with open(inp) as f:
		lines = f.readlines()
		n = 0
		es = []
		while n < len(lines):
			if 'END OF SELF-CONSISTENCY' in lines[n]:
				break
			n += 1
		while n < len(lines):
			if '------------------------' in lines[n]:
				n += 1
				break
			n += 1
		while n < len(lines):
			if '========================' in lines[n]:
				break
			e = float(lines[n][24:39])
			es.append(e)
			n += 1
		return es

def create_pp_and_basis(blochl_mult, sine_proj):
	psp = read_psp('PSP_SMOOTHED.DAT')
	fs = read_hfj_dft(os.path.join(BASIS_DIR, 'HFJ.DAT'))
	if os.path.exists(os.path.join(PROJ_DIR, 'HFJ.DAT')) and os.path.exists(os.path.join(PROJ_DIR, 'HFJ.RES')):
		projs = read_hfj_dft(os.path.join(PROJ_DIR, 'HFJ.DAT'))
		depth = get_hfj_depth(os.path.join(PROJ_DIR, 'hfj.inp'))
		es = get_hfj_energies(os.path.join(PROJ_DIR, 'HFJ.RES'))
		es = [e + depth for e in es]
		assert es and len(es) == len(projs), "different func count in HFJ.DAT and HFJ.RES!!!"
		projs = [p for p, e in zip(projs, es) if e > PROJ_E_LIMIT]
		psp.set_projectors(projs)

	psp.convert_potentials(True)
	sep_g = psp.create_separable_form(blochl_mult, sine_proj)
	psp.convert_potentials(False)
	sep_s = psp.create_separable_form(blochl_mult, sine_proj)

	basis = Basis()
	basis.load_ze_from_pp(sep_g)
	basis.set_functions(fs)
	basis.ortonorm()

	return sep_g, sep_s, basis

