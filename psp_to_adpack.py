#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.util.log_colorizer import init_log
from qcldm.applications.psp_sep import create_pp_and_basis
from qcldm.openmx_format.pao_format import PAO
from qcldm.openmx_format.vps_format import VPS

init_log(sys.argv)

blochl_mult = int(sys.argv[1])
sine_num = int(sys.argv[2])

sepg, seps, basis = create_pp_and_basis(blochl_mult, sine_num)

name = sepg.name()
VPS.to_file(sepg, '%s_GRECP.vps' % name)
VPS.to_file(seps, '%s_SL.vps' % name)
PAO.to_file(basis, '%sG%.1f.pao' % (name, basis.cutoff))








