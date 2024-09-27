#!/usr/bin/python
import re, sys, math, logging
import itertools
sys.dont_write_bytecode = True

from qcldm.turbomole_format.turbo_basis import TurboBasis
from qcldm.atom.shells import Shells


template = '''\\begin{table*}[!h]
\\small
\\caption{Some ECPs}
\\label{table:generated_ecps} 
\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}%s}
\\hline
\\hline
%s
\\end{tabular*}
\\end{table*}'''

ecp_names =  ['%s %s' % (e,n) for e, n in zip(*(iter(sys.argv[1:]),) * 2)]
#print(ecp_names)
turbo_basises, turbo_ecps = TurboBasis.read_basis('basis')
body = '\\multicolumn{%d}{c}{ECPs, ECPs, ECPs...} \\\\\n\\hline\n' % (len(ecp_names) * 4)

for en in ecp_names:
	assert en in turbo_ecps.keys(), 'No ecp with name \'%s\' in basis file' % en

body += ' & '.join(['\multicolumn{4}{c%s}{%s}' % ('' if (i == len(ecp_names) - 1) else '|', n.replace('_', '\\_')) for i, n in enumerate(ecp_names)]) + '\\\\\n'
tmp = []
for j, n in enumerate(ecp_names):
	tmp.extend(['\multicolumn{4}{c%s}{n$\mathrm{_{core}}$ = %d; $\mathrm{l_{max}}$ = %d}' % ('' if (j == len(ecp_names) - 1) else '|', turbo_ecps[n].ncore, len(turbo_ecps[n].semilocal))])
body += ' & '.join(tmp) + '\\\\\n'
body += ' & '.join(['&  \multicolumn{1}{c}{Exponent}', '\multicolumn{1}{c}{r$\\mathrm{^n}$}', '\multicolumn{1}{c|}{Coefficient}'] * (len(ecp_names) - 1) + ['&  \multicolumn{1}{c}{Exponent}', '\multicolumn{1}{c}{r$\\mathrm{^n}$}', '\multicolumn{1}{c}{Coefficient}']) + '\\\\\n\\hline\n'

lmax = max([len(turbo_ecps[x].semilocal) for x in ecp_names])
table = []
for i in range(lmax + 1):
	tc = []
	for j in range(len(ecp_names) * 4):
		tc.append([])
	table.append(tc)
for i in range(lmax + 1):
	for j, n in enumerate(ecp_names):
		ecp = turbo_ecps[n]
		ecp_lmax = len(ecp.semilocal)
		ls = ' '
		part = None
		if i == 0:
			ls = Shells.SHELLS[ecp_lmax].lower()
			part = ecp.local
		elif i < len(ecp.semilocal) + 1:
			ls = '%s-%s' % (Shells.SHELLS[i - 1].lower(), Shells.SHELLS[ecp_lmax].lower())
			part = ecp.semilocal[i - 1]
		else:
			ls = ' '
		table[i][j * 4].append(ls)
		if part:
			for (k, f) in part.functions:
				table[i][j * 4 + 1].append('%.5f' % f.a)
				table[i][j * 4 + 2].append('%d' % f.l)
				table[i][j * 4 + 3].append('%.5f' % k)

for i, c in enumerate(table):
	body += ' & ' * (len(ecp_names) * 4 - 1) + '\\\\\n'
	data = itertools.zip_longest(*table[i], fillvalue='')
	for row in data:
		body += ' & '.join(row) + '\\\\\n'



with open('basis.tex', 'w') as outp:
	outp.write(template % ('|'.join(['lrrr'] * len(ecp_names)), body))










































	





