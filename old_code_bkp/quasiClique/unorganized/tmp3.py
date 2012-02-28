from networkx import *
from networkx.readwrite import *
from networkx.spectrum import *
from math import sqrt
import psyco
from scipy import sparse

psyco.full()

g = read_gpickle('yybPM8N7Q16R2.WUblast_cut30.step2ed.gpickle')
g.remove_all_selfloops()
g.remove_all_multiedges()
g.ban_multiedges()

n = g.number_of_nodes()
node_mapping = {}
node_mapping_count = 0
asp = sparse.lil_matrix((n,n))
for e in g.edges_iter():
	try:
		i = node_mapping[e[0]]
	except:
		node_mapping[e[0]] = node_mapping_count
		i = node_mapping_count
		node_mapping_count += 1
	try:
		j = node_mapping[e[1]]
	except:
		node_mapping[e[1]] = node_mapping_count
		j = node_mapping_count
		node_mapping_count += 1
	print "{0}\t{1}\t{2}".format(i,j,e[2][2])
#	asp[i,j] = e[2][2]
#	asp[j,i] = e[2][2]

