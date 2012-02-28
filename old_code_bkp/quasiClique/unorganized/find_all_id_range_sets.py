import os,re,sys,fnmatch
from networkx import XGraph
from networkx.readwrite.gpickle import *
from networkx.cliques import *
from Bio import SeqIO

def clique_find_next_member(members, candidates, G, desired_size):
	if len(candidates) == 0:
		return members

	for i in xrange(len(candidates)):
		c = candidates[i]
		members.append(c)
		finaled = clique_find_next_member(members,filter(lambda x: G.has_edge(c,x), candidates[(i+1):]),G,desired_size)
		if len(finaled) >= desired_size:
			return finaled
		members.pop()
	return members			

def find_max_pid_range_set(pickle_filename,range_min,range_max,max_size=15):
	"""
		Given a pid graph pickle, read it, then output largest set of seq IDs(nodes)
		s.t. the pairwise PID is between range <range_min>-<range_max>

		This can be done by first removing all edges that are < <range_min> 
  		OR > <range_max> then call the networkx clique functions
	"""
	G = read_gpickle(pickle_filename)
	edges_to_delete = filter(lambda e: e[2] < range_min or e[2] > range_max, G.edges_iter())
	G.delete_edges_from(edges_to_delete)
	print >> sys.stderr, "{0}: {1} nodes, {2} valid edges".format(pickle_filename,G.number_of_nodes(),G.number_of_edges())

	for node in G.nodes_iter():
		print >> sys.stderr, "find cliques containing node {0}".format(node)
		a_clique = clique_find_next_member([node], G.neighbors(node), G, max_size)
		if len(a_clique) >= max_size:
			return a_clique
	return []


dirname='/homes/gws/lachesis/Larry_Riboswitch/Latest_ribo_20081002_plus20071102_curated'
max_size = 15

for filename in fnmatch.filter(os.listdir(dirname),'*.sto'):
	best_clique = find_max_pid_range_set(dirname+'/pid_graphs/'+filename+'.pid_graph.pickle',.4,.6,max_size)
	rec_dict = SeqIO.to_dict(SeqIO.parse(open(dirname+'/fnas/'+filename[:-4]+'.fna'),'fasta'))	
	for id in best_clique[:max_size]:
		print ">{0}\n{1}".format(id,rec_dict[id].seq.tostring())
		sys.stdout.flush()

