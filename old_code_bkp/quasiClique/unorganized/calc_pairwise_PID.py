import os,re,sys,fnmatch
from miscAlignment import pairwise_pid
from Bio import SeqIO
from networkx import XGraph
from networkx.readwrite import *

"""
	Calculate pairwise percent identity for all pairs of ncRNAs for each family,
	then storing them as a networkx XGraph, where the edge weight is PID.
"""

dirname = '/homes/gws/lachesis/Larry_Riboswitch/Latest_ribo_20081002_plus20071102_curated'
for filename in fnmatch.filter(os.listdir(dirname),'*.sto'):
	print >> sys.stderr, "{0}......".format(filename)
	with open(os.path.join(dirname,filename)) as handle:
		recs = list(SeqIO.parse(handle, "stockholm"))

	pid_graph = XGraph()
	for i,r1 in enumerate(recs):
		for r2 in recs[(i+1):]:
			pid_graph.add_edge(r1.id,r2.id,pairwise_pid(r1.seq.tostring(),r2.seq.tostring()))

	write_gpickle(pid_graph, os.path.join(dirname, filename+'.pid_graph.pickle'))
