import math, sys
import matplotlib.pyplot as plt
import numpy as np


def dfac(n):
	return 1 if n < 2 else reduce(lambda x,y: y*x, list(range(n,1,-2)))

def gauss_norm(a, l):
	return (2**(2*l+3.5)  / dfac(2*l+1) / math.pi**0.5)**0.5 * a**((2.*l+3)/4)


class GaussFunction:
	def __init__(self, a, l):
		self.a = a
		self.l = l

	def __call__(self, r):
		return (r**self.l) * math.exp(-self.a * r**2)

class GaussFunctionContracted:
	def __init__(self):
		self.fs = []

	def __call__(self, r):
		res = 0
		for c, f in self.fs:
			res += f(r) * c
		return res

class EcpPart:
	def __init__(self):
		self.base = None
		self.add = None

def read_ecp(fname, number):
	with open(fname) as f:
		lines = list(f)
		n = 0
		while n < len(lines) - 1:
			ls = lines[n].split()
			if len(ls) == 2 and ls[0] == str(number) and lines[n+1].strip() == 'INPUT':
#				print(n, lines[n].strip())
				break
			n += 1
		if n == len(lines) -1:
			return
		l_sizes = [int(x) for x in lines[n+2].split()[1:]]
		max_l = -1
		for l_size in l_sizes:
			if l_size != 0:
				max_l += 1
			else:
				break
#		print(l_sizes)
		n = n + 3
		parts = {}
		for l in range(len(l_sizes)):
			if l_sizes[l] == 0:
				continue
			part = EcpPart()
			part.base = GaussFunctionContracted()
			part.add = GaussFunctionContracted()
			
#			print('================')
			for line in lines[n:n+l_sizes[l]]:
				ls = line.split()
#				print(line.strip())
				ga = float(ls[0])
				gc = float(ls[1])
				gl = int(ls[2]) + 2
				g = GaussFunction(ga, gl)
				if len(ls) == 3:
					part.base.fs.append([gc, g])
				elif len(ls) == 4:
					part.add.fs.append([gc, g])
			real_l = max_l if l == 0 else l - 1
			parts[real_l] = part
			n += l_sizes[l]
		
		for l in parts.keys():
			if l == max_l:
				continue
			part = parts[l]
			part.base.fs.extend(parts[max_l].base.fs)
			part.add.fs.extend(parts[max_l].add.fs)
		return parts
			
		
				
parts = {}			
parts[0] = read_ecp(sys.argv[1] + '.start', int(sys.argv[2]))
parts[1] = read_ecp(sys.argv[1] + '.restart', int(sys.argv[2]))
assert(parts[0].keys() == parts[1].keys())

fig, axs = plt.subplots(len(parts[1]), 2)

for n in range(2):
	for l in sorted((list(parts[1].keys()))):
		x = np.logspace(-10, 1, num=500)
		axs[l, n].set_xscale('log')
		axs[l, n].plot(x, [parts[n][l].base(xx) for xx in x], color='blue', linewidth=0.5)
		axs[l, n].plot(x, [parts[n][l].add(xx) for xx in x], color='green', linewidth=0.5)
		axs[l, n].plot(x, [parts[n][l].base(xx) + parts[n][l].add(xx) for xx in x], color='red', linewidth=0.5)
		axs[l, n].axhline(y=0, color='black', linestyle='--', linewidth=0.5)
		axs[l, n].set_title('spdfgh'[l], loc='left')
plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1)	
fig.suptitle(sys.argv[1])
plt.show()

































