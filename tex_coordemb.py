#!/usr/bin/python
import re, sys, math, logging, os
import itertools
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.atom.shells import Shells
from qcldm.structures.atom_vector import AtomKeys

template = '''\\begin{table*}[!h]
\\small
\\caption{Cluster structure}
\\label{table:generated_embedding} 
\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}%s}
\\hline
\\hline
%s
\\end{tabular*}
\\end{table*}'''

paths = sys.argv[1:]
cfs = [ControlFormat.from_path(x) for x in paths]
body = ''
#body = '\\multicolumn{%d}{c}{Embedding} \\\\\n\\hline\n' % 5
body += ' & '.join(['\multicolumn{5}{c%s}{%s}' % ('' if (i == len(paths) - 1) else '|', n.replace('_', '\\_')) for i,n in enumerate(paths)]) + '\\\\\n'
body += ' & '.join(['\multicolumn{1}{c}{Atom}', '\multicolumn{1}{c}{x \AA}', '\multicolumn{1}{c}{y \AA}', '\multicolumn{1}{c}{z \AA}', '\multicolumn{1}{c|}{partial charge}'] * (len(cfs) - 1) + ['\multicolumn{1}{c}{Atom}', '\multicolumn{1}{c}{x \AA}', '\multicolumn{1}{c}{y \AA}', '\multicolumn{1}{c}{z \AA}', '\multicolumn{1}{c}{partial charge}']) + '\\\\\n\\hline\n'
body += ' & '.join([''] * 10) +  '\\\\\n'
table = []

for j in range(len(cfs) * 5):
	table.append([])
for i, cf in enumerate(cfs):
	for a in cf.cell.atoms:
		table[5 * i].append(a.name())
		table[5 * i + 1].append('%.5f' % a.position().x)
		table[5 * i + 2].append('%.5f' % a.position().y)
		table[5 * i + 3].append('%.5f' % a.position().z)
		table[5 * i + 4].append('%.5f' % a.data()[AtomKeys.ESTIMATED_CHARGE] if AtomKeys.ESTIMATED_CHARGE in a.data().keys() else '--')

data = itertools.zip_longest(*table, fillvalue='')
for row in data:
	body += ' & '.join(row) + '\\\\\n'


with open('coord.tex', 'w') as outp:
	outp.write(template % ('|'.join(['lrrrr'] * len(cfs)), body))





































	





