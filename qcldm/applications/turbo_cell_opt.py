#!/usr/bin/python 
import re, os, sys, threading, time
from math3d import Vector
import numpy as np
import fortranformat as ff
from ..structures.atom_vector import AtomVector
from ..crystal_format.crystal_out import CrystalOut
from ..crystal_format.crystal_part_writer import write_crystal_part
from ..turbomole_format.turbo_grad_utils import read_atoms, write_atoms, read_grad, write_grad, clear_relax
from ..util.xyz_format import write_xyz
from functools import reduce

class OptData:
	def __init__(self):
		cell_file = 'cell.start'
		coord_file = 'coord.start'
		if os.path.exists('cell.restart') and os.path.exists('coord.restart'):
			print('restarting')
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

	def cleared_x(self, x):
		res = []
		for i, a in enumerate(self.atoms):
			if a.name() == 'zz':
				res.extend([0.,0.,0.])
			else:
				res.extend(x[3*i:3*i+3])
		return res

	def coord_delta(self, old_atoms, new_atoms, damp):
		assert len(old_atoms) == len(new_atoms)
		res = []
		for o, n in zip(old_atoms, new_atoms):
			res.append((n.position().x - o.position().x) * damp)
			res.append((n.position().y - o.position().y) * damp)
			res.append((n.position().z - o.position().z) * damp)
		return res

	def apply_delta(self, dx):
		for i, a in enumerate(self.atoms):
			v = Vector(np.array(dx[i*3:i*3+3]) + a.position()._data)
			a.set_position(v)

	def dscf_step(self):
		write_atoms(self.atoms, 'coord')
		assert os.system('part_dscf > log_dscf.log') == 0
		time.sleep(1)

	def grad_step(self):
		assert os.system('part_grad > log_grad.log') == 0
		time.sleep(1)

	def read_grads(self):
		x = read_grad(len(self.atoms))
		xr = self.real_x(x)
		xv1, _ = self.co.cell.cut_and_expand_by_vectors(xr, self.real_atoms(), self.atoms)
		xc1, _ = self.co.cell.cut_and_expand_by_coords(xr, self.real_atoms(), self.atoms)
		xv1r = self.real_x(xv1)
		xc1r = self.real_x(xc1)
		grad = (reduce(lambda s, i: s + i**2, xr, 0.))**0.5
		grad_v = (reduce(lambda s, i: s + i**2, xv1r, 0.))**0.5
		grad_c = (reduce(lambda s, i: s + i**2, xc1r, 0.))**0.5	
		return grad, grad_v, grad_c	

	def rewrite_grad(self, vector_step):
		x = read_grad(len(self.atoms))
		xr = self.real_x(x)
		x1, _ = self.co.cell.cut_and_expand_by_vectors(xr, self.real_atoms(), self.atoms) if vector_step else self.co.cell.cut_and_expand_by_coords(xr, self.real_atoms(), self.atoms)
		write_grad(self.cleared_x(x1))

	def rewrite_coords(self, vector_step, damp):
		old_atoms = self.atoms
		new_atoms = read_atoms('coord')
		dx = self.coord_delta(old_atoms, new_atoms, damp)
		dxr = self.real_x(dx)
		dx1, _ = self.co.cell.cut_and_expand_by_vectors(dxr, self.real_atoms(), self.atoms) if vector_step else self.co.cell.cut_and_expand_by_coords(dxr, self.real_atoms(), self.atoms)
		self.co.cell.modify_by_vectors(dxr, self.real_atoms()) if vector_step else self.co.cell.modify_by_coords(dxr, self.real_atoms())
		self.apply_delta(dx1)

	def relax_step(self, damp, vector_step):
		self.rewrite_grad(vector_step)
		assert os.system('part_relax > log_relax.log') == 0
		time.sleep(1)
		self.rewrite_coords(vector_step, damp)
		write_crystal_part(self.co, 'cell.restart')
		write_atoms(self.atoms, 'coord.restart')
		
	def write_result(self, grad, grad_v, grad_c, vector_step, damp):
		typ = 'VCTR' if vector_step else 'ICRD'
		with open("coord.log", "a") as log:
			log.write("%s | grad = %11.8f | grad_v = %11.8f | grad_c = %11.8f | damp = %4.2f\n" % (typ, grad, grad_v, grad_c, damp))

	def write_result_init(self, grad, grad_v, grad_c, damp):
		typ = 'INIT'
		with open("coord.log", "a") as log:
			log.write("%s | grad = %11.8f | grad_v = %11.8f | grad_c = %11.8f | damp = %4.2f\n" % (typ, grad, grad_v, grad_c, damp))
		
	def semistep(self, damp, vector_step):
		if self.first:
			self.dscf_step()
			self.grad_step()
			grad, grad_v, grad_c = self.read_grads()
			self.write_result_init(grad, grad_v, grad_c, damp)
			self.first = False
		self.relax_step(damp, vector_step)
		self.dscf_step()
		self.grad_step()
		grad, grad_v, grad_c = self.read_grads()
		self.write_result(grad, grad_v, grad_c, vector_step, damp)
		return grad_v if vector_step else grad_c

	def step(self, damp):
		grad_c = self.semistep(damp, False)
		clear_relax()
		grad_v = self.semistep(damp, True)
		clear_relax()
		return max(grad_c, grad_v)
		
	def optimize_alternating(self, eps):
		grad = 999999999999
		damp = 0.1
		while True:
			grad1 = self.step(damp)
			if abs(grad1-grad) < eps:
				break
			grad = grad1
		
	def optimize_longrun(self, eps):
		grad = 999999999999
		damp = 0.1
		while True:
			while True:
				grad1 = self.semistep(damp, False)
				if abs(grad1-grad) < eps:
					break
				grad = grad1
			clear_relax()
			while True:
				grad1 = self.semistep(damp, True)
				if abs(grad1-grad) < eps:
					break
				grad = grad1
			clear_relax()
			if grad < eps * 10:
				break

	def optimize(self, eps):
		with open("coord.log", "w") as log:
			log.write("OPTIMIZATION START: eps=%.2e\n" % eps)
#		self.optimize_alternating(eps)
		self.optimize_longrun(eps)
		with open("coord.log", "a") as log:
			log.write("OPTIMIZATION FINISH\n")






























