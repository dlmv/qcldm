import os, re

from ..util.fileutils import make_dir
from ..record_format.psp_processing import smooth_psp, psp_to_hfj
from ..record_format.psp_format import read_psp
from ..functions.numeric_function import NumericOperations
from ..atom.shells import Shells
from ..util.elements import ELEMENTS

hfj_inp_template = '''
 %NAME%
 KL =   0
 NS =  %BASIS_COUNT%
 NSO=   0
 DFT=   2
 KEC= 101 130
 SIC=   1
 KAV=   0
 KPR=   0
 Z  = %ZVAL%
 LAG=   1
 JM =  -2.00
 RWL=   %CUTOFF% %DEPTH%
 R2 =   %GRIDEND%

         NL   J       QQ     KP   QQ_frac

%INPUT%


========================================
'''[1:-1]

BASIS_DIR = "BASIS"
PROJ_DIR = "PROJECTORS"

MAXN = 15
MAXL = 3

PROJ_N = 5

VIRTUAL_OCC = 0.001
CHARGE = 1.

LINE_FORMAT = " %3d    %2d%s (%d/%d)   %6.4f    %d    %5.3f\n"

def prepare_psp(limit, cutoff, depth = -20):
	func = lambda f: NumericOperations.parabolic_smooth(f, limit)
	smooth_psp('PSP.DAT', 'PSP_SMOOTHED.DAT', func)

	make_dir(BASIS_DIR)
	psp_to_hfj('PSP_SMOOTHED.DAT', os.path.join(BASIS_DIR, 'HFJ.POT'))

	make_dir(PROJ_DIR)
	psp_to_hfj('PSP_SMOOTHED.DAT', os.path.join(PROJ_DIR, 'HFJ.POT'), True)

	psp = read_psp('PSP_SMOOTHED.DAT')

	t = hfj_inp_template
	t = t.replace("%NAME%", psp.name())
	t = t.replace("%ZVAL%", "%.2f" % (-psp.zval))

	i = 1
	inplist = ""
	eld = ELEMENTS[psp.z].eleconfig_dict
	nval = 0
	for n, l in list(eld.keys()):
		if n > nval:
			nval = n
	valence = Shells.estimate_valence(psp.z)
	for n in range(1, MAXN + 1):
		for l in range(MAXL + 1):
			n_min = min(psp.pps[l].keys())
			real_n = n_min + n - 1
			occ = VIRTUAL_OCC * (4. * l + 2)
			key = (real_n, Shells.SHELLS[l].lower())
			if key in eld:
				occ = eld[key]
				val = real_n + max(0, l - 1) == nval
				if val:
					occ *= (valence - CHARGE) / valence
			for jj in range(abs(2*l-1), 2*l+2, 2):
				jocc = occ / (4. * l + 2) * (jj + 1)
				line = LINE_FORMAT % (i, n + l, Shells.SHELLS[l], jj, 2, jocc, 0, 1)
				inplist += line
				i += 1

	t = t.replace("%INPUT%", inplist)


	with open(os.path.join(BASIS_DIR, 'hfj.inp'), 'w') as f:
		tb = t
		tb = tb.replace("%BASIS_COUNT%", str(MAXN * (2 * MAXL + 1)))
		tb = tb.replace("%CUTOFF%", "%.2f" % (cutoff))
		tb = tb.replace("%DEPTH%", "%.1f" % (-abs(depth)))
		tb = tb.replace("%GRIDEND%", "%.6f" % (cutoff - 1e-6))
		f.write(tb)

	with open(os.path.join(PROJ_DIR, 'hfj.inp'), 'w') as f:
		tb = t
		tb = tb.replace("%BASIS_COUNT%", str(PROJ_N * (2 * MAXL + 1)))
		tb = re.sub("\s*RWL.*", "", tb)
#		tb = tb.replace("%CUTOFF%", "%.2f" % (psp.grid[-1] * 10))
#		tb = tb.replace("%DEPTH%", "%.1f" % (-0))
		tb = tb.replace("%GRIDEND%", "%.6f" % (-psp.grid[-1]))
		f.write(tb)

























	
