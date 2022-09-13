import re, sys, math, os, logging
from ..functions.numeric_function import NumericFunction, NumericOperations
from .record_format import read_records, write_records
from .psp_format import read_psp_headers, read_pp

#smoothing

def write_pp(data, zval):
	return [(v * r**2 + zval * r) for r, v in data]

def smooth_psp(inpname, outpname, smooth_func):
	records = read_records(inpname)
	npf = records[0][0][1]
	ngrid = int(records[0][0][2])
	raw_pps = read_psp_headers(records[0][0], npf)
	grid = records[1][0][:ngrid]
	for i in range(int(npf)):
		if raw_pps[i][1] >= 0:
			data, _, zval = read_pp(records[i + 2], grid)
			data = smooth_func(data)
			records[i + 2][0][:len(data)] = write_pp(data, zval)
		else:
			pass
	write_records(records, outpname)

#psp to hfj

def rewrite_headers(raw_pps, data):
	nl = {}
	raw_pps = [p for p, k in raw_pps if k >= 0]
	for p in raw_pps:
		if p.l not in list(nl.keys()):
			nl[p.l] = (p.n, p.n)
		else:
			nl[p.l] = (min(nl[p.l][0], p.n), max(nl[p.l][1], p.n))

	for i, p in enumerate(raw_pps):
		pos = 20 + i * 6
		data[pos] = p.n - nl[p.l][0] + p.l + 1
		data[pos + 1] = p.l
		data[pos + 2] = nl[p.l][1] - p.n
		data[pos + 3] = 0
		data[pos + 4] = 0
		data[pos + 5] = p.j * 2

	return nl

def load_pp(record, grid, pp, pf):
	record[0] = list(record[0])
	record[1] = list(record[1])
	pp.data = list(zip(grid, record[0][:len(grid)]))
	pf.data = list(zip(grid, record[1][:len(grid)]))
	pp.zval = float(record[0][len(grid) + 10])
	pf.gam = record[1][len(grid) + 1]
	pp.gam = record[0][len(grid) + 1]
	pf.sp = record[0][len(grid)]
	pf.coef = record[1][len(grid) + 3:len(grid) + 7] + record[1][len(grid) + 8:len(grid) + 9]
	pp.coef = record[0][len(grid) + 4:len(grid) + 8]

def psp_to_hfj(inpname, outpname, localonly=False):
	records = read_records(inpname)
	newrecords = []
	header = [0] * len(records[0][0])
	header[0:20] = records[0][0][0:20]

	z = records[0][0][0]

	npot = int(records[0][0][1])
	header[18] = npot * 2 + 2

	ngrid = int(records[0][0][2])
	raw_pps = read_psp_headers(records[0][0], npot)
	grid = records[1][0][:ngrid]

	nl = rewrite_headers(raw_pps, header)

	newrecords = [[header, header]]

	newrecords.append(records[1])

	pps = []
	pfs = []
	for i in range(int(npot)):
		if raw_pps[i][1] >= 0:
			p = raw_pps[i][0]
			f = NumericFunction()
			f.n = p.n
			f.l = p.l
			f.j = p.j
			load_pp(records[i + 2], grid, p, f)

			pps.append(p)
			pfs.append(f)

	pps_map = NumericOperations.resort_lnj(pps)
	lmax = sorted(pps_map.keys())[-1]
	nmax_lmax = sorted(pps_map[lmax].keys())[-1]
	loc = NumericOperations.j_to_ae(pps_map[lmax][nmax_lmax][lmax - 0.5], pps_map[lmax][nmax_lmax][lmax + 0.5])[0]
	loc.coef = [((c2 * lmax + c1 * (lmax + 1)) / (2. * lmax + 1)) for c2, c1 in zip(pps_map[lmax][nmax_lmax][lmax - 0.5].coef, pps_map[lmax][nmax_lmax][lmax + 0.5].coef)]

	new_pps = []
	for p in pps:
		np = p
		if p.l != lmax or p.n != nmax_lmax:
			np = p - loc
			np.coef = [c1 - c2 for c1, c2 in zip(p.coef, loc.coef)]
			if localonly:
				np *= 0
				np.coef = [0] * len(p.coef)
		np.zval = p.zval
		new_pps.append(np)

	pps = new_pps
	new_pps = []
	pps_map = NumericOperations.resort_lnj(pps)

	for p in pps:
		np = p
		if p.n != nl[p.l][1]:
			np = p - pps_map[p.l][nl[p.l][1]][p.j]
			np.coef = [c1 - c2 for c1, c2 in zip(p.coef, pps_map[p.l][nl[p.l][1]][p.j].coef)]
		np.zval = p.zval
		new_pps.append(np)

	pps = new_pps
	new_pps = []
	pps_map = NumericOperations.resort_lnj(pps)

	for p in pps:
		np = p
		if p.l != 0 and p.n == nl[p.l][1]:
			if p.j == p.l - 0.5:
				np = NumericOperations.j_to_ae(pps_map[p.l][p.n][p.l - 0.5], pps_map[p.l][p.n][p.l + 0.5])[0]
				np.coef = [((c2 * p.l + c1 * (p.l + 1)) / (2. * p.l + 1)) for c2, c1 in zip(pps_map[p.l][p.n][p.l - 0.5].coef, pps_map[p.l][p.n][p.l + 0.5].coef)]
			elif p.j == p.l + 0.5:
				np = NumericOperations.j_to_ae(pps_map[p.l][p.n][p.l - 0.5], pps_map[p.l][p.n][p.l + 0.5])[1]
				np.coef = [(c1 - c2) for c2, c1 in zip(pps_map[p.l][p.n][p.l - 0.5].coef, pps_map[p.l][p.n][p.l + 0.5].coef)]
			else:
				assert False
		np.j = p.j
		np.zval = p.zval
		new_pps.append(np)

	for p in new_pps:
		record = [0] * len(records[0][0])
		record[0:ngrid] = [v for r, v in p.data]
		for i in reversed(list(range(len(p.data)))):
			if p.data[i][1] != 0:
				record[len(p.data) + 2] = i + 1
				break
		record[len(p.data) + 3] = -2
		record[len(p.data) + 4:len(p.data) + 8] = p.coef
		record[len(p.data) + 8] = z
		record[len(p.data) + 10] = p.zval
		newrecords.append([record, record])

	for f in pfs:
		record = [0] * len(records[0][0])
		record[0:ngrid] = [v for r, v in f.data]
		for i in reversed(list(range(len(f.data)))):
			if f.data[i][1] != 0:
				record[len(f.data) + 2] = i + 1
				break
		record[len(f.data)] = f.sp
		record[len(f.data) + 1] = 1.0
		record[len(f.data) + 3] = f.gam
		record[len(f.data) + 4:len(f.data) + 9] = f.coef
		newrecords.append([record, record])

	write_records(newrecords, outpname)




