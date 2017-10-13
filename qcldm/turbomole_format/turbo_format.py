import re, os, sys

class TurboTemplate:

	class Param:
		def __init__(self, key, lineparam, multiparam):
			self.key = key
			self.lineparam = lineparam
			self.multiparam = multiparam

	def __init__(self, data):
		self.data = data

	key_regex = re.compile(u"\$[\w\-]+")

	def param(self, key):
		for d in self.data:
			if d.key == key:
				return d
		return None

	@staticmethod
	def from_string(datastring):
		lines = datastring.splitlines()
		data = []
		buf = []
		key = ''
		lineparam = ''
		for l in lines + ["$end"]:
			m = TurboTemplate.key_regex.match(l)
			if m:
				if key:
					p = TurboTemplate.Param(key, lineparam, buf)
					data.append(p)
				key = m.group()[1:]
				lineparam = l[len(key) + 1:]
				buf = []
			elif key:
				buf.append(l)
			else:
				print "Strange line: ", l

		f = TurboTemplate(data)
		return f

	@staticmethod
	def from_file(name):
		with open(name) as f:
			return TurboTemplate.from_string(f.read())

	def to_string(self):
		res = ''
		for d in self.data:
			res += "$%s%s\n" % (d.key, d.lineparam)
			for l in d.multiparam:
				res += "%s\n" % l
		return res

	def to_file(self, name):
		with open(name, 'w') as f:
			f.write(self.to_string())

