#!/usr/bin/python

class AtomKeys:
	BADER_CHARGE = 'Bader Charge'
	MULLIKEN_CHARGE = 'Mulliken Charge'
	ORBITAL_COUNT = 'orb_count'
	CUTOFF = 'cutoff'
	FULL_VALENCE = 'full_valence'
	ESTIMATED_VALENCE = 'estimated_valence'
	ORBITAL_ARRAY = 'orbital_array'
	ESTIMATED_CHARGE = 'estimated_charge'

class AtomVector:

	def __init__(self, name, vec, data=None):
		self._name = name
		self._vector = vec
		self._data = data or {}

	def position(self):
		return self._vector

	def name(self):
		return self._name

	def data(self):
		return self._data

	def relative(self, center):
		return AtomVector(self._name, self._vector - center._vector, self._data)

	def distance(self, other):
		return self.position().dist(other.position())
