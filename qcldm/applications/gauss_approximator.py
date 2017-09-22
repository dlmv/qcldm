import re, sys, math, logging
sys.dont_write_bytecode = True

from scipy.optimize import lsq_linear, fmin_slsqp
from ..functions.gauss_function import GaussFunctionContracted, GaussFunctionNormed
import copy

MAX_REC = 2
GRID_K = 4 #grid size is 2k+1

def expand_grid(f):
	grid = f.grid()
	m = grid[-1] / grid[-2]
	d = grid[-1] - grid[-2]
	cutoff = grid[-1]
	while grid[-1] < cutoff * 25:
		grid.append(grid[-1] * m)
	return grid

def approximate(f, conv, logbound, maxcoef, nmin, nmax, writefile=False):
	logging.info(u'')
	logging.info(u'*********************************************')
	logging.info(u'  Approximating a function, l=%d, n=%d' % (f.l, f.n))
	logging.info(u'*********************************************')
	logging.info(u'')
	fnorm = f.norm()
	logging.info(u'Norm = %g' % (fnorm))
	f.normalize()
	nn = f.num_nodes()
	logging.info(u'Num nodes = %d' % (nn))
	f = f.rescaled(expand_grid(f))
	convd = False
	minconv = (-1, 1e50, None)

	for n in range(nmin, nmax + 1):
		logging.info(u'Assume n = %d, Conv = %e' % (n, conv))
		cg, convd, res = approximate_n(f, n, logbound, conv, maxcoef, writefile)
		fg = cg.to_numeric(f.grid())

		#Wrong norm - skip
		if (cg.norm() > fg.norm() * 1.05):
			convd = False
			continue

		#Wrong zero behaviour - skip
		if abs(cg(f.data[0][0]) - f.data[0][1]) > 0.01:
			convd = False
			continue

		#Also skip
		for r, v in f.data:
			if abs(cg(r) - v) > 1:
				convd = False
				continue

		#First good step, continue
		if not minconv[2]:
			minconv = (n, res, cg)
			convd = False
			continue

		#Converged, but fast descent - continue
		if convd and (res <= minconv[1] / 1.05):
			convd = False

		#Check if better
		if res < minconv[1] and ((abs(cg.norm() - 1) < abs(minconv[2].norm() - 1)) or abs(cg.norm() - 1) < res * 15):
			minconv = (n, res, cg)

		if convd:
			break
	n, res, cg = minconv
	logging.info(u'Taking n = %d, Conv = %e' % (n, res))
	
	fg = cg.to_numeric(f.grid())

	logging.info(u'AnalyticNorm = %f' % (cg.norm()))
	logging.info(u'NumericNorm = %f' % (fg.norm()))

	cg.fs = [(k*fnorm, g) for k, g in cg.fs]

	logging.info(u'RestoredNorm = %f' % (cg.norm()))

#	cg, res = shift_exps(f, cg, GRID_K, conv, maxcoef, res)

#	logging.info(u'New conv = %e' % (res))

#	fg = cg.to_numeric(f.grid())

#	logging.info(u'AnalyticNorm = %f' % (cg.norm()))
#	logging.info(u'NumericNorm = %f' % (fg.norm()))

	if writefile:
		print_approx('approx_l%d_n%d' % (f.l, f.n), f, cg)
	return cg, res

def approximate_n(f, n, logbound, conv, maxcoef, writefile):
	fun = lambda x: approx_res(f, x[0], x[1], n, conv, maxcoef)
	res, loga, logarn = scan_2grid(fun, -logbound / 2, logbound / 2, logbound / 2, GRID_K, conv)
	a = math.exp(loga)
	r = (math.exp(logarn) / a)**(1. / n)
	cg = approx_func(f, a, r, n, conv, maxcoef)
	convd = res < conv
	logging.info(u'%s: %e, Norm = %f' % ("Converged" if convd else "Not converged", res, cg.norm()))
	if writefile:
		print_approx('approx_tmp', f, cg)
	return cg, convd, res


def scan_2grid(fun, xs1, xs2, r, k, conv, recurdepth=1):
	mins = []
	ratio = 2
	res0 = fun( (xs1, xs2) )
	indent = ' ' * recurdepth
	logging.debug(indent + u'Depth %d. Start: %f (%f, %f)' % (recurdepth, res0, xs1, xs2))
	for i1 in range(-k + 1, k):
		for i2 in range(-k + 1, k):
			x1 = xs1 + r * i1 / k
			x2 = xs2 + r * i2 / k
			res = fun( (x1, x2) )
			if not mins:
				mins = [(res, x1, x2)]
			else:
				ins = True
				for (mres, mx1, mx2) in mins:
					if mres * ratio < res:
						ins = False
						break
				if ins:
					mins.append((res, x1, x2))
					mins = [item for item in mins if item[0] < res * 2]
					mins.sort(key = lambda x: x[0])
					if len(mins) > k:
						mins = mins[:k]
	res, x1, x2 = mins[0]
	logging.debug(indent + u'Current result: %f (%f, %f)' % (res, x1, x2))
	if res > max(res0 * 0.95, ((res0 + conv) / 2)):
		logging.debug(indent + u'No need to search here!')
		return res, x1, x2
	recurdepth += 1
	if recurdepth > MAX_REC:
		logging.debug(indent + u'Recursion limit!')
		return res, x1, x2
	optmins = []
	for (res, x1, x2) in mins:
		optmins.append(scan_2grid(fun, x1, x2, r / k, k, conv, recurdepth))
		if optmins[-1][0] < conv:
			break
	optmins.sort(key = lambda x: x[0])
	return optmins[0]



def gauss_set(l, a, r, n):
	progression = [a * r**i for i in range(n)]
	gs = []
	for p in progression:
		gs.append(GaussFunctionNormed(p, l))
	return gs

def gauss_approx(f, gs, conv, maxcoef):
	b = []
	a = []
	for x, v in f.data:
		b.append(v)
		ta = []
		for g in gs:
			ta.append(g(x))
		a.append(ta)
	return lsq_linear(a, b, (-maxcoef, maxcoef), 'bvls', conv/100, None, None, None, 0)

#def shift_exps(f, cg, n, conv, maxcoef, res0):
#	logging.info(u'Tweaking result...')
#	logstep = 1.
#	gs = [g for k, g in cg.fs]
#	rmin = res0
#	while True:
#		rstart = rmin
#		for i in range(len(gs)):
#			logstart = math.log(gs[i].a) - logstep if i == 0 else math.log(gs[i-1].a)
#			logend = math.log(gs[i].a) + logstep if i == len(gs) - 1 else math.log(gs[i+1].a)
#			step = (logend - logstart) / (2 * n)
#			for k in range(1, 2 * n):
#				newa = math.exp(logstart + k * step)
#				newg = GaussFunctionNormed(newa, gs[i].l)
#				newgs = gs[:i] + [newg] + gs[i + 1:]
#				res = gauss_approx(f, newgs, conv, maxcoef)
#				r = (sum(map(lambda x:x*x,res.fun)) / len(res.fun))**0.5

#				tcg = GaussFunctionContracted()
#				for g, x in zip(gs, res.x):
#					tcg.fs.append((x, g))
#				fg = tcg.to_numeric(f.grid())

#				#Wrong norm - skip
#				if (tcg.norm() > fg.norm() * 1.05) or abs(tcg.norm() - 1) > 0.3:
#					continue

#				#Wrong zero behaviour - skip
#				if abs(tcg(f.data[0][0]) - f.data[0][1]) > 0.01:
#					continue

#				#Also skip
#				for r, v in f.data:
#					if abs(tcg(r) - v) > 1:
#						convd = False
#						continue

#				#Check if better
#				if r < rmin and ((abs(tcg.norm() - 1) < abs(cg.norm() - 1)) or abs(tcg.norm() - 1) < r * 20):
#					cg = tcg
#					gs = newgs
#					rmin = r

#		logging.debug(u'Current result: %e' % rmin)
#		if rmin > rstart * 0.95:
#			logging.debug(u'No need to search here!')
#			break
#	return cg, rmin

#def f_tmp(f, gs, ns):
#	cg = GaussFunctionContracted()
#	cg.fs = zip(ns, gs)
#	res = 0
#	for (r, fv) in f.data:
#		gv = cg(r)
#		res += (fv - gv)**2
#	return res

#def norm_tmp(gs, ns):
#	cg = GaussFunctionContracted()
#	cg.fs = zip(ns, gs)
#	return (cg.norm() - 1)**2

#def gauss_approx_slsqp(f, gs, conv, maxcoef):
#	cg = GaussFunctionContracted()
#	for g in gs:
#		cg.fs.append((1., g)) 
#	cg.normalize()
#	initial = [n for n, g in cg.fs]
#	func = lambda ns: f_tmp(f, gs, ns)
#	normfunc = lambda ns: norm_tmp(gs, ns)
#	bounds = [(-maxcoef, maxcoef)] * len(gs)
#	res = fmin_slsqp(func, initial, [normfunc], None, [], None, [], None, None, None, (), 1000, 1e-10, 0, None, 1)
#	print '================'
#	print res
#	return res


def approx_res(f, loga, logarn, n, conv, maxcoef):
	a = math.exp(loga)
	r = (math.exp(logarn) / a)**(1. / n)
	gs = gauss_set(f.l, a, r, n)
	res = gauss_approx(f, gs, conv, maxcoef)
	return (sum(map(lambda x:x*x,res.fun)) / len(res.fun))**0.5
#	res = gauss_approx_slsqp(f, gs, conv, maxcoef)
#	print (res[1]/ len(res[0]))**0.5
#	return (res[1]/ len(res[0]))**0.5

def approx_func(f, a, r, n, conv, maxcoef):
	gs = gauss_set(f.l, a, r, n)
	res = gauss_approx(f, gs, conv, maxcoef)
	cg = GaussFunctionContracted()
	for g, x in zip(gs, res.x):
		cg.fs.append((x, g))
	return cg
#	res = gauss_approx_slsqp(f, gs, conv, maxcoef)
#	cg = GaussFunctionContracted()
#	for g, x in zip(gs, res.fx):
#		cg.fs.append((x, g))
#	return cg

def print_approx(name, f, cg):
	fg = cg.to_numeric(f.grid())
	outp = open(name, 'w')
	line = "%15.10f %15.10f %15.10f\n"
	for (x, v), (x1, v1) in zip(f.data, fg.data):
		outp.write(line % (x, v, v1))

def get_res(f, cg):
	res = 0
	fg = cg.to_numeric(f.grid())
	for (x, v), (x1, v1) in zip(f.data, fg.data):
		res += (v - v1)**2
	return res



