

def frange(x, y, jump):
	while x < y:
		yield x
		x += jump

def ufu(u, f):
	return u.getH().dot(f).dot(u)
