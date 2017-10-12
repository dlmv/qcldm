import re, os, sys


class Openmx_format:

	class RawText:

		def __init__(self, strings):
			self.strings = strings

	class Param:

		def __init__(self, key, value):
			self.key = key
			self.value = value

	class MultiParam:

		def __init__(self, key, value):
			self.key = key
			self.value = value

	def __init__(self, data):
		self.data = data

	def set_param(self, key, value):
		for d in self.data:
			if isinstance(d, Openmx_format.Param):
				if d.key == key:
					d.value = value
					return
		p = Openmx_format.Param(key, value)
		self.data.append(p)

	def set_multiparam(self, key, value):
		for d in self.data:
			if isinstance(d, Openmx_format.MultiParam):
				if d.key == key:
					d.value = value
					return
		p = Openmx_format.MultiParam(key, value)
		self.data.append(p)

	def param(self, key):
		for d in self.data:
			if isinstance(d, Openmx_format.Param):
				if d.key == key:
					return d.value
		return None

	def multiparam(self, key):
		for d in self.data:
			if isinstance(d, Openmx_format.MultiParam):
				if d.key == key:
					return d.value
		return None

	def check_param(self, d, keys):
		if isinstance(d, Openmx_format.Param) or isinstance(d, Openmx_format.MultiParam):
			for k in keys:
				if d.key == k:
					return True
		return False

	def remove_params(self, keys):
		check = lambda x: not self.check_param(x, keys)
		self.data[:] = [x for x in self.data if check(x)]

	def retain_params(self, keys):
		check = lambda x: self.check_param(x, keys)
		self.data[:] = [x for x in self.data if check(x)]

	def add_text(self, text):
		self.data.append(Openmx_format.RawText([text]))
		#TODO: add text before/after param

	@staticmethod
	def from_string(datastring):
		lines = datastring.splitlines()
		data = []
		buf = []
		key = ''
		for l in lines:
			if l.find('#') != -1:
				l = l[:l.find('#')]
			ls = filter(None, re.split("\s+", l.strip()))
			if not ls:
				if not data or not isinstance(data[-1], Openmx_format.RawText):
					data.append(Openmx_format.RawText(['']))
			elif key:
				if len(ls) == 1 and ls[0][-1] == '>':
					p = Openmx_format.MultiParam(key, buf)
					data.append(p)
					key = ''
				else:
				 	buf.append(ls)
			elif len(ls) == 1 and ls[0][0] == '<' and '.' in ls[0]:
					key = ls[0][1:]
					buf = []
			elif '.' in ls[0] or ls[0] == 'AtomSpecies':
				p = Openmx_format.Param(ls[0], ls[1:])
				data.append(p)
			else:
				if not data or not isinstance(data[-1], Openmx_format.RawText):
					data.append(Openmx_format.RawText(['']))

		f = Openmx_format(data)
		return f

	@staticmethod
	def from_file(name):
		with open(name) as f:
			return Openmx_format.from_string(f.read())

	def to_string(self):
		res = ''
		for d in self.data:
			if isinstance(d, Openmx_format.Param):
				res += ('{:35}' + ' {}'*len(d.value) + '\n').format(*tuple([d.key] + d.value))
			elif isinstance(d, Openmx_format.MultiParam):
				res += '<{}\n'.format(d.key)
				ml = 0
				for ls in d.value:
					ml = max(ml, len(ls))
				lengths = [0] * ml
				for ls in d.value:
					for i in range(len(ls)):
						if len(ls) > i:
							l = ls[i]
							if l[0] != '-':
								l = '-' + l
							lengths[i] = max(lengths[i], len(l))
				for ls in d.value:
					tmp = ''
					for i in range(len(ls)):
						tmp += ('{:>' + str(lengths[i]) + '} ').format(ls[i])
					tmp += "\n"
					res += tmp
				res += '{}>\n'.format(d.key)
			elif isinstance(d, Openmx_format.RawText):
				for s in d.strings:
					res += '{}\n'.format(s)
		return res

	def to_file(self, name):
		with open(name, 'w') as f:
			f.write(self.to_string())


