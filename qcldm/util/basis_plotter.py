import sys
sys.dont_write_bytecode = True
import math
from PIL import Image
import colorsys



point = [(0,0), (0,-1), (0,1), (-1,0), (1,0)]

def plot_funcs(fs, filename, xmax=-1, ymax = 2):
	h = 800
	w = 1200
	img = Image.new('RGB', (w, h), '#aaaaaa')
	pix = img.load()

	xmin = 0
	if xmax < 0:
		xmax = round(fs[0].data[-1][0])
	ymin = -ymax
#	for f in fs:
#		for r, v in f.data:
#			if v > ymax:
#				ymax = v
#			if v < ymin:
#				ymin = v
	imgdata = list(img.getdata())
	y0 = int(round(-ymin  * h / (ymax - ymin)))
	if y0 >= 0 and y0 <= h:
		for x in range(w):
			pix[x, y0] = (0,0,0)
		for np in range(int(xmax)):
			for yt in range(-h / 25, h / 25):
				pix[round(1.0 * (np - xmin) / (xmax - xmin) * w), y0 + yt] = (0,0,0)
	for nf in range(len(fs)):
		f = list(reversed(fs))[nf]
		hue = 5.0 / 6 * (len(fs) - nf - 1) / len(fs)
		color = tuple([int(round(i * 255)) for i in colorsys.hsv_to_rgb(hue, 1, 1)])

		xp = None
		yp = None
		for r, v in f.data:
			x = round(1.0 * (r - xmin) / (xmax - xmin) * w)
			y = round(1.0 * (-v - ymin) / (ymax - ymin) * h)
			if xp != None:
				n = int(round(max(abs(x - xp), abs(y - yp))))
				if n > 0 and n < 1000:
					for i in xrange(n + 1):
						xc = int(round((x * i + xp * (n - i)) / n))
						yc = int(round((y * i + yp * (n - i)) / n))
						for xs, ys in point:
							xx = xc + xs
							yy = yc + ys
							if xx >= 0 and xx < w and yy >= 0 and yy < h:
								pix[xx, yy] = color
			xp = x
			yp = y

	img.save(filename + '.png')

def cut_hue(v):
	pass

def plothue():
	h = 200
	w = 1200
	img = Image.new('RGB', (w, h), '#aaaaaa')
	pix = img.load()
	for x in range(w):
		hue = 1.0 * x / w
		color = tuple([int(round(i * 255)) for i in colorsys.hls_to_rgb(hue, 0.5, 1)])
		for y in range(h) if color.count(0) == 2 else range(h / 2):
			pix[x, y] = color
	img.save('hue.png')


#plothue()











