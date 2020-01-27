#!/usr/bin/python
import re, sys, math, logging
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
write_xyz(co.cell.supercell, 'supercell.xyz')

write_crystal_part(co, 'tmp')
co = CrystalOut.from_file('tmp')
write_crystal_part(co, 'tmp1')


