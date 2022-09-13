import re, sys, math, os, logging

from .record_format import read_records
from ..functions.numeric_function import NumericFunction, NumericOperations
from ..atom.pseudo_potential import GRECPPseudoPotential

def read_grid(data, n):
	return data[:n]

def read_pp(record, grid):
	pp = []
	pf = []
	zval = float(record[0][len(grid) + 10])
	pp = [( r, p/r**2 - zval / r) for r, p in zip(grid, record[0][:len(grid)])]
	pf = [( r, f / r) for r, f in zip(grid, record[1][:len(grid)])]
	return pp, pf, zval

def read_psp_headers(data, n):
	raw_pps = []
	for i in range(int(n)):
		p = NumericFunction()
		pos = 20 + i * 6
		p.n = int(data[pos])
		p.l = int(data[pos + 1])
		kp = data[pos + 3]
		p.j = data[pos + 5] / 2
		raw_pps.append((p, kp))
	return raw_pps

def read_psp(name):
	logging.info('')
	logging.info('*********************************************')
	logging.info('  Reading PSP file: %s' % name)
	logging.info('*********************************************')
	logging.info('')
	records = read_records(name)
	z = int(records[0][0][0])
	npf = records[0][0][1]
	ngrid = int(records[0][0][2])
	raw_pps = read_psp_headers(records[0][0], npf)
	grid = records[1][0][:ngrid]
	pps = []
	pfs = []
	zval = None
#	rcmap = {}
	for i in range(int(npf)):
		if raw_pps[i][1] >= 0:
			p = raw_pps[i][0]
			f = NumericFunction()
			f.n = p.n
			f.l = p.l
			f.j = p.j
			p.data, f.data, zvalp = read_pp(records[i + 2], grid)

			if zval != None:
				assert zval == zvalp, "Different zval in PSP file"
			zval = zvalp

			pps.append(p)
			pfs.append(f)
		else:
			pass
	logging.debug('Number of PPs: %d' % len(pps))
	for p in pps:
		logging.debug(' L=%d N=%d J=%3.1f' % (p.l, p.n, p.j))
	pps_map = NumericOperations.resort_lnj(pps)
	lmax = sorted(pps_map.keys())[-1]
	nmax_lmax = sorted(pps_map[lmax].keys())[-1]
	logging.debug('Local PP: L=%d N=%d' % (lmax, nmax_lmax))
	loc = NumericOperations.j_to_ae(pps_map[lmax][nmax_lmax][lmax - 0.5], pps_map[lmax][nmax_lmax][lmax + 0.5])[0]

	gp = GRECPPseudoPotential(pps, pfs, loc, z, zval)
	return gp

















