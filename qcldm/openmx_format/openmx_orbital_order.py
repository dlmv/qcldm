from ..atom.harmonics import Orbital_Order

class Openmx_Order(Orbital_Order):

	def order(self, l):
		if l == 1:
			return [2,0,1]
		elif l == 2:
			return [2,4,0,3,1]
		elif l == 3:
			return [3,4,2,5,1,6,0]
		else:
			return Orbital_Order.order(self,l)


