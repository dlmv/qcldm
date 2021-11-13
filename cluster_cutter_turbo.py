#!/usr/bin/python
import re, sys, math, logging
sys.dont_write_bytecode = True

from qcldm.turbomole_format.control_format import ControlFormat
from qcldm.util.log_colorizer import init_log
from qcldm.util.xyz_format import write_xyz
from qcldm.embedding.cluster import Cluster
from qcldm.embedding.embedding_settings import EmbeddingSettings
from qcldm.structures.bader_reader import read_baders
from qcldm.structures.atom_vector import AtomKeys

init_log(sys.argv)



settingsfile = sys.argv[1]
settings = EmbeddingSettings.from_file(settingsfile)

c = ControlFormat().from_file('control', settings.reverse_embedding_map)
#for a in c.cell.atoms:

		

cluster = Cluster(c.cell, settings)
cluster.estimate_charges_dumb()
cluster.write()

#c = CrystalMeta()
#c.load('.')

#co = CrystalOut.from_file(c.out_file)
#write_xyz(co.cell.cell, 'cell.xyz')
#write_xyz(co.cell.supercell, 'supercell.xyz')

#try:
#	read_baders(co.cell)
#except Exception:
#	pass

#settingsfile = sys.argv[1]
#settings = EmbeddingSettings.from_file(settingsfile)


#cluster = Cluster(co.cell, settings)
#cluster.estimate_charges_dumb()
#cluster.write()





