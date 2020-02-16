#!/usr/bin/python 
import re, os, sys, threading, time
from math3d import Vector
import numpy as np
import fortranformat as ff
from ..structures.atom_vector import AtomVector
from ..crystal_format.crystal_out import CrystalOut
from ..crystal_format.crystal_part_writer import write_crystal_part
from ..turbomole_format.turbo_grad_utils import read_atoms, write_atoms, read_grad, write_grad
from ..util.xyz_format import write_xyz

class OptData:
	def __init__(self):
		cell_file = 'cell.start'
		coord_file = 'coord.start'
		if os.path.exists('cell.restart') and os.path.exists('coord.restart'):
			print 'restarting'
			cell_file = 'cell.restart'
			coord_file = 'coord.restart'
		self.co = CrystalOut.from_file(cell_file)
		self.atoms = read_atoms(coord_file)
		self.first = True
		write_xyz(self.co.cell.supercell, 'supercell.xyz')
		
	def real_atoms(self):
		res = []
		for a in self.atoms:
			if a.name() != 'zz':
				res.append(a)
		return res
		
	def real_x(self, x):
		res = []
		for i, a in enumerate(self.atoms):
			if a.name() != 'zz':
				res.extend(x[3*i:3*i+3])
		return res

	def coord_delta(self, old_atoms, new_atoms):
		assert len(old_atoms) == len(new_atoms)
		res = []
		for o, n in zip(old_atoms, new_atoms):
			res.append(n.position().x - o.position().x)
			res.append(n.position().y - o.position().y)
			res.append(n.position().z - o.position().z)
		return res

	def apply_delta(self, dx, damp):
		for i, a in enumerate(self.atoms):
			v = Vector(np.array(dx[i*3:i*3+3]) * damp + a.position()._data)
			a.set_position(v)

	def dscf_step(self):
		write_atoms(self.atoms, 'coord')
		assert os.system('/home/demidov/turbo/bin/amd64/part_dscf > log_dscf.log') == 0
		time.sleep(1)

	def grad_step(self, vector_step):
		assert os.system('/home/demidov/turbo/bin/amd64/part_grad > log_grad.log') == 0
		time.sleep(1)
		x = read_grad(len(self.atoms))
		xr = self.real_x(x)
		x1, _ = self.co.cell.cut_and_expand_by_vectors(xr, self.real_atoms(), self.atoms) if vector_step else self.co.cell.cut_and_expand_by_coords(xr, self.real_atoms(), self.atoms)
		write_grad(x1)
		x1r = self.real_x(x1)
		grad = (reduce(lambda s, i: s + i**2, x1r, 0.))**0.5
		return grad

	def relax_step(self, damp, vector_step):
		old_atoms = self.atoms
		assert os.system('/home/demidov/turbo/bin/amd64/part_relax > log_relax.log') == 0
		time.sleep(1)
		new_atoms = read_atoms('coord')
		dx = self.coord_delta(old_atoms, new_atoms)
		dxr = self.real_x(dx)
		dx1, _ = self.co.cell.cut_and_expand_by_vectors(dxr, self.real_atoms(), self.atoms) if vector_step else self.co.cell.cut_and_expand_by_coords(dxr, self.real_atoms(), self.atoms)
		self.co.cell.modify_by_vectors(dxr, self.real_atoms()) if vector_step else self.co.cell.modify_by_coords(dxr, self.real_atoms())
		self.apply_delta(dx1, damp)
		write_crystal_part(self.co, 'cell.restart')
		write_atoms(self.atoms, 'coord.restart')
		
	def semistep(self, damp, vector_step):
		if self.first:
			self.dscf_step()
			self.grad_step(vector_step)
			self.first = False
		self.relax_step(damp, vector_step)
		self.dscf_step()
		grad = self.grad_step(vector_step)
		return grad
		
	def step(self, damp):
		grad_c = self.semistep(damp, False)
		self.write_result(grad_c, 'coord', damp)
		grad_v = self.semistep(damp, True)
		self.write_result(grad_v, 'vector', damp)
		return max(grad_c, grad_v)
		
	def optimize(self, eps):
		grad = None
		with open("coord.log", "w") as log:
			log.write("OPTIMIZATION START: eps=%.2e\n" % eps)
		damp = 0.5
		while True:
			grad1 = self.step(damp)
			if not grad:
				grad = grad1
				continue
			if abs(grad1-grad) < eps:
				break
			grad = grad1


	def write_result(self, grad, typ, damp):
		with open("coord.log", "a") as log:
			log.write("grad = %13.8f | step = %s | damp = %5.3f\n" % (grad, typ, damp))




























