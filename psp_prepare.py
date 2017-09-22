#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.util.log_colorizer import init_log
from qcldm.applications.psp_preparer import prepare_psp

init_log(sys.argv)

cutoff = int(sys.argv[1])
limit = int(sys.argv[2])

prepare_psp(limit, cutoff)










