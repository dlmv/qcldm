import re, sys, math, os, logging

from record_format import read_records
from ..functions.numeric_function import NumericFunction

def read_common_header(records):
	z = records[0][0][0]
	nvf = int(records[0][0][1])
	ngrid = int(records[0][0][2])
	grid = records[1][0][:ngrid]
	return grid, nvf

def read_function(data, grid):
	return  [( r, f / r) for r, f in zip(grid, data[:len(grid)])]

#====PSF========

def read_psf_headers(data, n):
	funcs = []
	for i in range(int(n)):
		f = NumericFunction()
		pos = 20 + i * 6
		f.n = int(data[pos])
		f.l = int(data[pos + 1])
#		jvf = data[pos + 2]
#		res = data[pos + 3]
#		nodes = data[pos + 4]
		f.j = data[pos + 5] / 2
		funcs.append(f)
	return funcs

def read_psf(name):
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Reading PSF file: %s' % name)
	logging.info(u'*********************************************')
	logging.info(u'')
	records = read_records(name)
	grid, nvf = read_common_header(records)
	funcs = read_psf_headers(records[0][0], nvf)
	for i in range(int(nvf)):
		funcs[i].data = read_function(records[i + 2][0], grid)
	logging.debug(u'Number of functions: %d' % len(funcs))
#	for p in funcs:
#		logging.debug(u' L=%d N=%d J=%3.1f' % (p.l, p.n, p.j))
	return funcs

#=======HFJ-DFT=========

def read_hfj_dft_headers(data, n):
	funcs = []
	for i in range(int(n)):
		f = NumericFunction()
		pos = 20 + i * 3
		f.n = int(data[pos])
		kk = int(data[pos + 1])
#		qq = data[pos + 2]
		f.l = kk if kk > 0 else -kk - 1
		f.j = 1.*abs(kk) - 0.5
		funcs.append(f)
	return funcs

def read_hfj_dft(name):
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Reading HFJ-DFT output: %s' % name)
	logging.info(u'*********************************************')
	logging.info(u'')
	records = read_records(name)
	grid, nvf = read_common_header(records)
	funcs = read_hfj_dft_headers(records[0][0], nvf)
	for i in range(nvf):
		funcs[i].data = read_function(records[i + 4][0], grid)
	logging.debug(u'Number of functions: %d' % len(funcs))
#	for p in funcs:
#		logging.debug(u' L=%d N=%d J=%3.1f' % (p.l, p.n, p.j))
	return funcs

#======HFD=======

def read_hfd_headers(data, n):
	funcs = []
	for i in range(int(n)):
		f = NumericFunction()
		pos = 20 + i * 6
		f.n = int(data[pos])
		f.l = int(data[pos + 1])
#		qq = data[pos + 2]
#		qw = data[pos + 3]
#		kk = data[pos + 4][
		f.j = data[pos + 5] / 2
		funcs.append(f)
	return funcs

def read_hfd(name):
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Reading HFD output: %s' % name)
	logging.info(u'*********************************************')
	logging.info(u'')
	records = read_records(name)
	grid, nvf = read_common_header(records)
	funcs = read_hfd_headers(records[0][0], nvf)
	for i in range(nvf):
		funcs[i].data = read_function(records[i + 4][0], grid)
	logging.debug(u'Number of functions: %d' % len(funcs))
#	for p in funcs:
#		logging.debug(u' L=%d N=%d J=%3.1f' % (p.l, p.n, p.j))
	return funcs







