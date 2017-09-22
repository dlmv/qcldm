import re, sys, math, os, logging

import fortranformat as ff
import numpy as np

RECORD_SIZE = 512

def read_records(name):
	if name.lower().endswith('.dat') or name.lower().endswith('.pot'):
		return read_records_binary(name)
	elif name.lower().endswith('.num'):
		return read_records_text(name)
	else:
		raise RuntimeError(name)

def write_records(records, name):
	if name.lower().endswith('.dat') or name.lower().endswith('.pot'):
		write_records_binary(records, name)
	elif name.lower().endswith('.num'):
		write_records_text(records, name)
	else:
		raise RuntimeError(name)


def read_record_text(lines, n):
	record = [[], []]
	if n >= len(lines) or 'RECORD' not in lines[n]:
		return None, n
	n += 1
	while n < len(lines) and 'RECORD' not in lines[n]:
		ls = filter(None, re.split('\s*', lines[n]))
		record[0].append(float(ls[1]))
		record[1].append(float(ls[2]))
		n += 1
	return record, n

def read_records_text(name):
	f = open(name)
	lines = f.readlines()
	f.close()
	records = []
	n = 0
	while n < len(lines):
		r, n = read_record_text(lines, n)
		records.append(r)
	return records

def read_records_binary(name):
	records = np.fromfile(name)
	assert len(records) % (RECORD_SIZE*2) == 0, 'non integer number of records'
	records = [[records[x:x+RECORD_SIZE], records[x+RECORD_SIZE:x+RECORD_SIZE*2]] for x in xrange(0, len(records), RECORD_SIZE*2)]
	return records

def write_records_binary(records, name):
	records = [records[i][j][k] for i in xrange(0, len(records)) for j in xrange(2) for k in xrange(RECORD_SIZE)]
	np.array(records).tofile(name)

def write_records_text(records, name):
	f = open(name, 'w')
	for i, r in enumerate(records):
		f.write("     RECORD =  %2d\n" % (i + 1))
		for j, (a, b)  in enumerate(zip(r[0], r[1])):
			s = ff.FortranRecordWriter("( 1X,I11,1X,G24.16E3,1X,G24.16E3 )").write([j + 1, a, b])
			f.write(s + "\n")








