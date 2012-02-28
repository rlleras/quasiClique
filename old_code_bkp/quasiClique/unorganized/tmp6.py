import os,re,sys
from networkx import *
from networkx.readwrite import *
import cluster_steps as c
import pClique as p

def run(blast_filename):
	filename = c.step1_process_blast('',blast_filename,30)
	return
	g = read_gpickle(filename)
	g.delete_nodes_from(filter(lambda n: g.degree(n)<=2, g.nodes_iter()))
	g_nodes = g.nodes()
	print >> sys.stderr, "number of nodes...",g.number_of_nodes()
	S, H = p.convert_graph_connectivity_to_sparse(g, g_nodes)
	valid_rows = range(g.number_of_nodes())
	QQ = []
	last_QQ_size = 0
	failed_count = 0
	while failed_count < 10:
		p.get_cliques(S, H, QQ, 0.8, 10, valid_rows)
		if len(QQ) == last_QQ_size:
			failed_count += 1
		else:
			failed_count = 0
			last_QQ_size += 1
			
	f = open(blast_filename+'.gamma80iter10.cliques','w')
	for q in QQ:
		f.write(" ".join(map(lambda x: g_nodes[x], q))+"\n")
	print >> sys.stderr, "output written to....",f.name
	f.close()

if __name__=="__main__":
	d = sys.argv[1]
	for file in os.listdir(d):
		if not file.endswith('.gpickle') and not os.path.exists(d+'/'+file + '_cut30.step1ed.gpickle'):
			print >> sys.stderr, file
			run(d+'/'+file)
			
