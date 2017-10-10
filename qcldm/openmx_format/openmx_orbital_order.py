import numpy
from ..atom.quantum_numbers import Basis_Matrix

class Openmx_orbital_order(Basis_Matrix):

	def lmatrix(self, l):
		if l == 1:
			return numpy.matrix([[0,0,1],[1,0,0],[0,1,0]])
		elif l == 2:
			return numpy.matrix([[0,0,1,0,0],[0,0,0,0,1],[1,0,0,0,0],[0,0,0,1,0],[0,1,0,0,0]])
		elif l == 3:
			return numpy.matrix([[0,0,0,0,0,0,1],[0,0,0,0,1,0,0],[0,0,1,0,0,0,0],[1,0,0,0,0,0,0],[0,1,0,0,0,0,0],[0,0,0,1,0,0,0],[0,0,0,0,0,1,0]])
		else:
			return Basis_Matrix.lmatrix(self,l)


