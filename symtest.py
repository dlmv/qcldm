#!/usr/bin/python
import re, sys, math, logging, random
sys.dont_write_bytecode = True

from qcldm.crystal_format.crystal_out import CrystalOut
from qcldm.crystal_format.crystal_meta import CrystalMeta
from qcldm.crystal_format.crystal_part_writer import write_crystal_part
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.embedding.cluster import Cluster
from qcldm.embedding.embedding_settings import EmbeddingSettings
from qcldm.structures.bader_reader import read_baders
from qcldm.structures.atom_vector import AtomKeys

init_log(sys.argv)

c = CrystalMeta()
c.load('.')

co = CrystalOut.from_file(c.out_file)
write_xyz(co.cell.cell, 'cell.xyz')
write_xyz(co.cell.bordered_cell, 'cell_b.xyz')
write_xyz(co.cell.supercell, 'supercell.xyz')

co.cell.build_vector_parameters()

b = [random.random() for x in range(9)]

#print b
#print co.cell.cut_by_vectors_symmetry(b)

write_crystal_part(co, 'tmp')

#co.cell.set_coord_basis(bas)

#bas = list(co.cell.get_coord_basis())
#co.cell.set_coord_basis(bas)

#write_crystal_part(co, 'tmp1')
#bas = [x + random.random() for x in bas]
#co.cell.set_coord_basis(bas)
#write_crystal_part(co, 'tmp2')

#co = CrystalOut.from_file('tmp')
#write_crystal_part(co, 'tmp1')


