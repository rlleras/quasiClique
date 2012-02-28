from cluster_steps_settings import *
import cluster_steps2 as c2
from miscUtils import groupby
"""
=======================================================
Step 7 of Quasi-Clique Algorithm
=======================================================
Looks at a particular homology graph and resulting cliques
and ask questions like:
		(1) What can we say about the near misses?
		(2) What can we say about the false inclusions?
		...these two evaluate our choice of gamma
"""

class CliqueSummary:
	highest_free_id = 1
	def __init__(self, fam, Q, TPs, missed):
		self.id   = CliqueSummary.highest_free_id
		CliqueSummary.highest_free_id += 1
		self.fam  = fam
		self.Q    = Q
		self.size = len(Q)
		self.TPs  = TPs
		self.numTPs = len(TPs)
		self.FPs  = filter(lambda q: q not in TPs, Q)
		self.numFPs = self.size - self.numTPs
		self.missed = missed
		self.nummissed = len(missed)
		self.name = "{fam}_{id}_{size}".format(fam=self.fam, id=self.id, size=self.size)

		self.TP2TP = 0
		self.FP2TP = 0
		self.FP2FP = 0
		self.TP2miss = 0

	def __str__(self):
		return \
			"NAME   : {0}\n".format(self.name) + \
			"FAMILY : {0}\n".format(self.fam) + \
			"TP     : {0}\n".format(self.numTPs) + \
			"FP     : {0}\n".format(self.numFPs) + \
			"Missed : {0}\n".format(self.nummissed) + \
			"^gamma : {0:.2f}\n".format(self.__gamma__()) + \
			"TP2TP  : {0:.2f}\n".format(self.TP2TP*1. / (self.numTPs*(self.numTPs-1)/2)) + \
			"FP2TP  : {0:.2f}\n".format(self.FP2TP*1. / max(1,self.numFPs*self.numTPs)) + \
			"TP2miss: {0:.2f}\n".format(self.TP2miss*1. / max(1,(self.numTPs*self.nummissed)))
	
	def __gamma__(self):
		return (self.TP2TP + self.FP2TP + self.FP2FP)*1. / (self.size*(self.size-1)/2)

def label_sets_for_nodes(cursor):
	"""
	Returns <sets_for_nodes> id --> ncRNA family (can be None)
	"""
	c2.CONN_FUNC = CONN_FUNC
	nodes_label = {}
	cursor.execute("SELECT i FROM sets_for_nodes")
	for r in cursor.fetchall():
		print >> sys.stderr, "labeling for node....", r[0]
		ncRNA_id, ncRNA_fam = c2.check_hit(r[0])
		nodes_label[r[0]] = ncRNA_fam
	return nodes_label

def false_inclusions(clique_pickle, output_dir):
	"""
	For each clique, determine the "false inclusion" (FPs),
	and see how well connected they are to the TPs in the clique
	"""
	c2.CONN_FUNC = CONN_FUNC
	with CONN_FUNC() as cursor:
#		nodes_label = label_sets_for_nodes(cursor)
#		f = open('tmp.pickle','w')
#		dump(nodes_label, f)
#		f.close()
		with open('tmp.pickle') as handle:
			nodes_label = load(handle)

		with open(clique_pickle) as handle:
			QQQ = load(handle)
		f = open(os.path.join(output_dir, 'REPORT'), 'a')
		f.write("#---------------------------------------------------\n")
		f.write("CLIQUE PICKLE: {0}\n".format(clique_pickle))
		for Q in QQQ:
			# members: <sets_for_nodes>  id --> ncRNA family
			members = dict(map(lambda q: (q, nodes_label[q]), Q))
			tallied = [(k, len(g)) for k,g in groupby(members.values())]
			tallied.sort(key=itemgetter(1))
			fam = tallied[-1][0]
			if fam is None:
				continue
			print >> sys.stderr, "this is a {0} clique......".format(fam)
			TPs = filter(lambda q: members[q]==fam, Q)
			ss = ",".join( map(str, Q) )
			cursor.execute("SELECT i1,i2 FROM parsed WHERE i1 IN ({bunch}) AND i2 IN ({bunch})".format(bunch=ss))
			G = Graph()
			for i1,i2 in cursor.fetchall():
				G.add_edge(i1, i2)
			# find potential near misses
			missed = filter(lambda k: nodes_label[k]==fam and k not in Q, nodes_label.iterkeys())
			ss = ",".join( map(str, missed) )
			if len(ss) > 1:
				for q in Q:
					print >> sys.stderr, "running externals with ....", q
					cursor.execute("SELECT i1,i2 FROM parsed WHERE (i1={q} AND i2 IN ({bunch}))\
							OR (i2={q} AND i1 IN ({bunch}))".format(q=q, bunch=ss))
					for i1,i2 in cursor.fetchall():
						G.add_edge(i1, i2)
		
			CS = CliqueSummary(Q=Q, fam=fam, TPs=TPs, missed=missed)
			output_prefix = os.path.join(output_dir, CS.name)
			if len(Q) < 30:
				networkx_graph_to_graphviz(G, CS, output_prefix=output_prefix+'_caseALL', case=-1, recordCS=True)
			else:
				networkx_graph_to_graphviz(G, CS, output_prefix=output_prefix+'_case0', case=0, recordCS=True)
				for case in xrange(1,3):
					networkx_graph_to_graphviz(G, CS, output_prefix=output_prefix+'_case'+str(case), case=case, recordCS=False)
			f.write('@@@@@@@@@@@@@@@@@@@@@\n' + str(CS))			
		f.close()

def networkx_graph_to_graphviz(G, CS, output_prefix, case, recordCS):
	with open(output_prefix+'.gv', 'w') as f:
		f.write("""graph {name} {{\n
	  	    edge [ dir=none ];\n
		    node [ shape=circle, style=filled, fontsize=2.0, height=0.1, width=0.1, fixedsize=true ];\n""".format(name=CS.name))
		TP_i, FP_i, missed_i = 1, 1, 1
		map_i_to_node = {}
		for n in G.nodes_iter():
			if n in CS.TPs:
				f.write("""  TP{i} [ colorscheme=brbg8, color=7 ];\n""".format(i=TP_i))
				map_i_to_node[n] = "TP{i}".format(i=TP_i)
				TP_i += 1
			elif n in CS.FPs:
				f.write("""  FP{i} [ colorscheme=brbg8, color=2 ];\n""".format(i=FP_i))
				map_i_to_node[n] = "FP{i}".format(i=FP_i)
				FP_i += 1
			else:
				f.write("""  miss{i} [ color=gold1 ];\n""".format(i=missed_i))
				map_i_to_node[n] = "miss{i}".format(i=missed_i)
				missed_i += 1
	#	print("""  {{ rank=same; {TPs};}}""".format(TPs=" ".join( map(lambda x:map_i_to_node[x], TPs) )))
	#	print("""  {{ rank=same; {FPs};}}""".format(FPs=" ".join( map(lambda x:map_i_to_node[x], FPs) )))
		for n1,n2 in G.edges_iter():
			if n1 in CS.TPs and n2 in CS.TPs:
				if recordCS:
					CS.TP2TP += 1
				color = 9
				if case == -1 or case == 0:
					f.write("""  {n1} -- {n2} [ colorscheme=set312, color={color}];\n""".format(n1=map_i_to_node[n1],\
							n2=map_i_to_node[n2],color=color))
			elif n1 in CS.Q and n2 in CS.Q:
				if recordCS:
					if n1 in CS.TPs or n2 in CS.TPs:
						CS.FP2TP += 1
					else:
						CS.FP2FP += 1
				color = 6
				if case == -1 or case == 1:
					f.write("""  {n1} -- {n2} [ colorscheme=set312, color={color}];\n""".format(n1=map_i_to_node[n1],\
							n2=map_i_to_node[n2],color=color))
			else:
				if n1 in CS.FPs or n2 in CS.FPs:
					color = 7
				else:
					if recordCS:
						CS.TP2miss += 1
					color = 10
				if case == -1 or case == 2:
					f.write("""  {n1} -- {n2} [ colorscheme=set312, color={color}];\n""".format(n1=map_i_to_node[n1],\
							n2=map_i_to_node[n2],color=color))
		f.write("}\n")		
			
def graph_ncRNA_connectivity(clique_pickle=None):
	"""
	Looks at table:parsed, and for each ncRNA of family F,
	calculate the ratio of connected F members/Nones,
	and do this for both immediately connected nodes
	and path-2 nodes.
	"""
	c2.CONN_FUNC = CONN_FUNC
	ncRNA_map = defaultdict(lambda: []) # fam --> list of <i>s
	dum = 0
	with CONN_FUNC() as cursor:
		cursor.execute("SELECT i FROM sets_for_nodes")
		for r in cursor.fetchall():
			print >> sys.stderr, "i is....", r[0]
			id,fam = c2.check_hit(r[0])
			if id is not None:
				ncRNA_map[fam].append(id)
				dum += 1
			if dum > 100: break
		for fam,list_of_i in ncRNA_map.iteritems():
			get_connectivity(list_of_i, cursor)

def get_connectivity(list_of_i, cursor):
	G = Graph()
	list_of_i.sort()
	for k,i1 in enumerate(list_of_i):
		ss = ",".join(list_of_i[(k+1):])
		cursor.execute("SELECT i2 \
				FROM sets_for_nodes \
				WHERE i1={i1} \
				AND i2 in ({i2_bunch})".format(i1=i1, i2_bunch=ss))
		for r in cursor.fetchall():
			G.add_edge(i1, r[0])
		print(G.edges())
		raw_input()

if __name__ == "__main__":
	CONN_FUNC = get_conn_Firm
	c2.CONN_FUNC = get_conn_Firm
	false_inclusions('../../output/output_cliques/ALLFirm_RefSeq25_m30s0_M8N7Q16R2W3E2_cut35.quasi80_ALL.pickle', \
			'../../DRAW_DIR')
