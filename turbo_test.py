#!/usr/bin/python
import re, sys, os, logging
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.util.log_colorizer import init_log

init_log(sys.argv)

ControlFormat.from_file('control')






