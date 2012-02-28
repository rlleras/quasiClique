from cluster_steps_settings import *
import cluster_steps1 as c1
import pClique as p
from random import choice
from Bio import SeqIO
import ushuffle

NEW_MIN_CLIQUE_SIZE = 3

def new_cluster_pipe(rfam_fam, shuffle_ratio):
	assert type(shuffle_ratio) is int
	output_prefix = "Rfam_{fam}_shuffle{X}X".format(fam=rfam_fam, X=shuffle_ratio)
	fasta_filename = output_prefix+'.fna'
	blast_output = "{input}.M8N7Q16R2W3E2.WUblast".format(input=fasta_filename)

	report_f = open(output_prefix+'.report', 'w')
	if not os.path.exists(blast_output):
		dummy_id = 0
		nodes_to_index = {}
		with open(fasta_filename, 'w') as f:
			with get_conn_ncRNA() as cursor:
				cursor.execute("select id,seq from Rfam_fasta where rfam_fam='{fam}' order by id".format(fam=rfam_fam))
				for _id,seq in cursor.fetchall():
					id = "TP{0}_{1}".format(dummy_id, _id)
					f.write(">{id}\n{seq}\n".format(id=id, seq=seq))
					nodes_to_index[id] = dummy_id
					dummy_id += 1
					ushuffle.shuffle(seq, len(seq), 2)
					for x in xrange(shuffle_ratio):
						id = "FP{0}_{1}".format(dummy_id, _id)
						f.write(">{id}\n{seq}\n".format(id=id, seq=ushuffle.shuffle2()))
						nodes_to_index[id] = dummy_id
						dummy_id += 1
		start_t = time.time()				
		# now blast it
		os.system("xdformat -n -o {input} {input}".format(input=fasta_filename))
		os.system("blastn -d {input} -i {input} -M 8 -N -7 -Q 16 -R 2 -E 2 \
				-W 3 -mformat 2 -cpus 4 -o {output}".format(input=fasta_filename, output=blast_output))
		report_f.write("(1)  BLAST TIME: {0} sec\n".format(time.time()-start_t))

		# now parse the blast
		nodes_to_index = c1.NodesToIndex(nodes_to_index, -1)
		G = Graph()
		c1.step1_process_blast(blast_output=blast_output,\
				score_cutoff=35, nodes_to_index=nodes_to_index, G=G, program='WU')
		print >> sys.stderr, "Homology graph has {0} nodes, {1} edges....".format(\
				G.number_of_nodes(), G.number_of_edges())
		c1.export_to_db(G, nodes_to_index, 0, blast_output)
		# convert nodes_to_index into dict nodes_ind --> acc id
		nodes_to_index = dict( map(lambda (x,y):(y,x), nodes_to_index.d.items()) )
		with open(blast_output+'.nodes_to_index', 'w') as handle:
			for ind,id in nodes_to_index.iteritems():
				handle.write("{0}\t{1}\n".format(ind,id))

	# read back the .parsed and .sets_for_nodes files
	G = Graph()
	sets_for_nodes = {}
	nodes_to_index = {}
	with open(blast_output+'.parsed') as handle:
		for line in handle:
			raw = map(int, line.strip().split('\t'))
			G.add_edge(raw[0],raw[1])
	with open(blast_output+'.sets_for_nodes') as handle:
		for line in handle:
			raw = map(int, line.strip().split('\t'))
			sets_for_nodes[raw[0]] = {'nodes_ind':raw[1],'start':raw[2],'end':raw[3]}
	with open(blast_output+'.nodes_to_index') as handle:
		for line in handle:
			raw = line.strip().split()
			nodes_to_index[int(raw[0])] = raw[1]

	tmp = len(filter(lambda x: nodes_to_index[sets_for_nodes[x]['nodes_ind']].startswith('FP'), G.nodes_iter()))
	report_f.write("(2)  AFTER parsing BLAST, graph has {0} negative control nodes, {1} TP nodes\n".format(tmp, G.number_of_nodes()-tmp))

	# remove low deg (< 3) nodes
	x = filter(lambda n: G.degree(n)<NEW_MIN_CLIQUE_SIZE, G.nodes_iter())
	while len(x) > 0:
		G.delete_nodes_from(x)
		x = filter(lambda n: G.degree(n)<NEW_MIN_CLIQUE_SIZE, G.nodes_iter())
	tmp = len(filter(lambda x: nodes_to_index[sets_for_nodes[x]['nodes_ind']].startswith('FP'), G.nodes_iter()))
	report_f.write("(3)  AFTER recursively removing nodes of degree < 3, graph has {0} negative control nodes, {1} TP nodes\n".format(tmp, G.number_of_nodes()-tmp))
	report_f.write("----------------------------------------------------------------------------\n")
	report_f.write("OUT\tDIR\tCLIQUE_SIZE\tSCANNED_TP\tSCANNED_FP\tCM_time\n")
	report_f.write("----------------------------------------------------------------------------\n")

	# for now just brute force....go through node by node as seeds
	dummy_round = 0
	while G.number_of_nodes()>=NEW_MIN_CLIQUE_SIZE and G.number_of_edges()>=NEW_MIN_CLIQUE_SIZE:
		# find perfect max cliques with a random starting node
		G_nodes = G.nodes()
		S,H = p.convert_graph_connectivity_to_sparse(G, G_nodes)
		tQ = p.grasp(S, H, gamma=1.0, maxitr=20, given_starting_node=None)
		Q = map(lambda x: G_nodes[x], tQ)
		if len(Q) < NEW_MIN_CLIQUE_SIZE: # delete these nodes
			G.delete_nodes_from(Q)
			continue
		# PERFECT CLIQUE SANITY TESTING, DELETE LATER
		for x in Q:
			for y in Q:
				if x!=y:
					print >> sys.stderr, "testing....", x,y
					try:
						assert G.has_edge(x,y)
					except:
						return Q,G
		Q.sort()
		print >> sys.stderr, "clique is...", Q
		prefix = output_prefix + str(dummy_round) + '_size' + str(len(Q)) + '_'
		dummy_round += 1
		start_t = time.time()
		scan_dir,scan_result = run_cmfinder(Q, nodes_to_index, sets_for_nodes, prefix, os.path.abspath(fasta_filename))
		cm_time = time.time()-start_t
		if scan_result is not None:
			outf = open(os.path.basename(scan_dir)+'.gv','w')
			outf.write("""graph test{
			edge [ dir=none ];
			node [ style=filled, fontsize=2.0, height=0.1, width=0.1, fixedsize=true ];
			""")
			# draw this graph
			scanned = {'TP':0, 'FP':0}
			for n in G.nodes_iter():
				id = nodes_to_index[sets_for_nodes[n]['nodes_ind']] # id is something like TP1_NC_XXXX.... or FP10_NC_XXXX....
				_id = id[:id.find('_')]
				shape = 'circle' if id.startswith('TP') else 'box'
				if n in Q:
					outf.write("{0} [color=dodgerblue1, shape={1}];\n".format(n, shape))
				elif id in scan_result:
					outf.write("{0} [color=darkorange, shape={1}];\n".format(n, shape))
					scanned[_id[:2]] += 1
				else:
					outf.write("{0} [color=grey, shape={1}];\n".format(n, shape))
			for (n1,n2) in G.edges_iter(data=False):
				id1 = nodes_to_index[sets_for_nodes[n1]['nodes_ind']]
				id1 = id1[:id1.find('_')]
				id2 = nodes_to_index[sets_for_nodes[n2]['nodes_ind']]
				id2 = id2[:id2.find('_')]
				outf.write("{0} -- {1};\n".format(n1,n2))
			outf.write("}")
			outf.close()
			report_f.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(outf.name,scan_dir,len(Q),scanned['TP'],scanned['FP'],cm_time))
			report_f.flush()
		# delete the edges from the graph
		G.delete_edges_from(itertools.combinations(Q, 2))
		# again, remove low-degree nodes
		x = filter(lambda n: G.degree(n)<NEW_MIN_CLIQUE_SIZE, G.nodes_iter())
		while len(x) > 0:
			G.delete_nodes_from(x)
			x = filter(lambda n: G.degree(n)<NEW_MIN_CLIQUE_SIZE, G.nodes_iter())  
	report_f.close()

def clique_perfects(G, seed_i):
	G_nodes = G.nodes()
	S,H = p.convert_graph_connectivity_to_sparse(G, G_nodes)
	bestQ = []
	for itr in xrange(PERFECT_MAXITR):
		alpha = uniform(.1, .9)
		Q = map(lambda x: G_nodes[x], p.construct(S, H, alpha, G_nodes.index(seed_i)))
		# NOTE: no guarantee that seed_i is in the clique
		if len(Q) > len(bestQ):
			bestQ = Q
	return bestQ

def run_cmfinder(Q, nodes_to_index, sets_for_nodes, prefix, fasta_filename):
	import tempfile, miscCMF
	old_dir = os.popen("pwd").read().strip()
	MOTIF_DIR = '../../MOTIF_DIR/experiment'
	prefix_dir = tempfile.mkdtemp(dir=MOTIF_DIR, prefix=prefix)
	
	F = FastaReader(fasta_filename)
	os.chdir(prefix_dir)
	handle = open( os.path.basename(prefix_dir)+'.fna', 'w' )
	for q in Q:
		m = sets_for_nodes[q]
		acc_id = nodes_to_index[m['nodes_ind']]
		handle.write(">{0}({1}-{2})\n".format(acc_id,m['start'],m['end']))
		handle.write("{0}\n".format(F[acc_id].seq[m['start']:(m['end']+1)]))
	handle.close()	

	# run CMfinder
	os.system("cmfinder.pl -def " + handle.name)
	# rank the motifs with rank_cmfinder.pl
	os.system("rank_cmfinder.pl -w -rank \"{0}.*.motif.*\" {0}.summary".format(handle.name))
	# run pfold_pscore, the summary is written to <fna_name>.pscore.summary
	# and in addition we get back a descending sorted list of (lod,motif)
	pscores = miscCMF.pfold_pscore(motif_dir='.')
	print >> sys.stderr, "prefix dir is", prefix_dir
	# take the highest pscore motif, cmbuild then cmsearch the original fasta
	# TODO: change behaviour later?
	scan_result = None
	if len(pscores) > 0:
		scan_result = run_infernal(motif=pscores[0][1], scan_filename=fasta_filename,\
				output_prefix=pscores[0][1]+'.Rfam')
		print >> sys.stderr, "scan result is...", scan_result
	print >> sys.stderr, "prefix dir is", prefix_dir
	os.chdir(old_dir)
	return os.path.basename(prefix_dir), scan_result

def run_infernal(motif, scan_filename, output_prefix):
	import miscInfernal
	# run infernal (NOTE THE HARD-CODED stuff!!! should change later)
	os.system("cmbuild {motif}.cmbuilded {motif}".format(motif=motif))
	os.system("cmsearch --fil-T-hmm 10 --toponly -T 10 " + \
			"--tabfile {out}.hits -o {out}.cmscanned ".format(out=output_prefix) + \
			"{motif}.cmbuilded {scan}".format(motif=motif, scan=scan_filename))

	result = miscInfernal.read_tab_cmscanned(tab_filename=output_prefix+'.hits')
	# for simplicity, remove from result hits that are motif members
	# currently the IDs are HARDCODED like:  >AARF01000040.1/2307-2550(1-244)
	for id in os.popen("grep \">\" {motif}.fasta".format(motif=motif)).read().strip().split('\n'):
		try:
			_id = id[1:id.find('(')]
			print >> sys.stderr, "removing from cmscan hits {0}....".format(_id)
			del result[_id]
		except:
			pass

	return result

if __name__=="__main__":
	crap = new_cluster_pipe(sys.argv[1], int(sys.argv[2]))

