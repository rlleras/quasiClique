import os,re,sys
sys.path.append('/home/etseng/UW_Larry_ncRNA/code/quasiClique/')

import cluster_steps1 as c1
from networkx import *

blast_split_filename = sys.argv[1]
output_i = int(sys.argv[2])

G = Graph()
nodes_to_index = c1.get_nodes_to_index('justPaeruginosa.IGRs_M8N7Q16R2W3S35E2.WU.nodes_to_index')
c1.step1_process_blast(blast_split_filename, 35, nodes_to_index, G)
c1.export_to_db(G, nodes_to_index, output_i, 'test')

