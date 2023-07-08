#!/usr/bin/python
import re, sys, math, logging
import itertools
sys.dont_write_bytecode = True

from qcldm.turbomole_format.turbo_basis import TurboBasis
from qcldm.atom.shells import Shells


template = '''\\begin{table*}[!h]
\\small
\\caption{Some basises}
\\label{table:generated_basises} 
\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}%s}
\\hline
%s
\\end{tabular*}
\\end{table*}'''

basis_names =  ['%s %s' % (e,n) for e, n in zip(*(iter(sys.argv[1:]),) * 2)]
#print(basis_names)
turbo_basises, turbo_ecps = TurboBasis.read_basis('basis')
body = '\\multicolumn{%d}{c}{basises, basises, basises...} \\\\\n\\hline\n' % (len(basis_names) * 2 + 1)

for bn in basis_names:
	assert bn in turbo_basises.keys(), 'No basis with name \'%s\' in basis file' % bn

body += '& ' + ' & '.join(['\multicolumn{2}{c%s}{%s}' % ('' if (i == len(basis_names) - 1) else '|', n.replace('_', '\\_')) for i,n in enumerate(basis_names)]) + '\\\\\n'
body += '& ' + ' & '.join(['\multicolumn{1}{c}{Exp}', '\multicolumn{1}{c|}{Coeff}'] * (len(basis_names) - 1) + ['\multicolumn{1}{c}{Exp}', '\multicolumn{1}{c}{Coeff}']) + '\\\\\n'


lmax = max([max([f.fs[0][1].l for f in turbo_basises[x].functions]) for x in basis_names])

table = []
linetable = []
for i in range(lmax + 1):
	tc = []
	for j in range(len(basis_names) * 2 + 1):
		tc.append([])
	table.append(tc)

for i in range(lmax + 1):
	ltc = []
	for j in range(len(basis_names)):
		ltc.append([])
	linetable.append(ltc)


for i in range(lmax + 1):
	ls = Shells.SHELLS[i].lower()
	table[i][0].append(ls)
	for j, n in enumerate(basis_names):
		basis = turbo_basises[n]
		part = [f for f in basis.functions if f.fs[0][1].l == i]
		for f in part:
			for k, g in f.fs:
				table[i][j * 2 + 1].append('%.5f' % g.a)
				table[i][j * 2 + 2].append('%.5f' % k)
				linetable[i][j].append(0)
			linetable[i][j][-1] = 1
		if linetable[i][j]:
			linetable[i][j][-1] = 0


for i, c in enumerate(table):
	body += '\\hline\n'
	body += ' & ' * (len(basis_names) * 2) + '\\\\\n'
	data = itertools.zip_longest(*table[i], fillvalue='')
	ldata = itertools.zip_longest(*linetable[i], fillvalue=0)
	for row, lrow in zip(data, ldata):
		body += ' & '.join(row) + '\\\\'
		for j,v in enumerate(lrow):
			if v:
				body += '\cline{%s-%s}' % (2*j+2, 2*j+3)
		body += '\n'
	if i != len(table) - 1:
		body += ' & ' * (len(basis_names) * 2) + '\\\\\n'



with open('basis.tex', 'w') as outp:
	outp.write(template % ('|'.join(['c'] + ['rr'] * len(basis_names)), body))


























	





