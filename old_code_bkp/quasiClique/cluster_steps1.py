from cluster_steps_settings import *

"""
===============================================================================
Step 1 of Quasi-Clique Algorithm
===============================================================================
Handles parsing of BLAST output into homology graphs and merging of overlapping 
entries in the database
"""

def update_db_node_sets(low, high):
	"""
	Take all from table:sets_for_nodes that have field:ind between <low>,<high>,
	merge any nodes that are from the same ind (see nodes_to_index) & overlap
	by at least 1 bp.

	For efficiency, splits the task in chunks of 1000.
	"""
	for i in xrange(low, high, +1000):
		update_db_node_sets_help(i, min(i+1000, high))
		
def update_db_node_sets_help(low, high):
	conn = CONN_FUNC()
	cursor = get_dict_cursor(conn)
	
	# first, get all the node index between <low>-<high>
	cursor.execute("SELECT ind "
			       "FROM nodes_to_index "
			       "WHERE {low} <= ind  and ind <= {high}".format(low=low, high=high))
	nodes_inds = map(lambda x: x['ind'], cursor.fetchall())
	print >> sys.stderr, "finished reading nodes_to_index........"

	# then examine, for each node index set, whether there are overlapping sets
	update_node_args = []
	update_parsed_args = []
	del_node_args = []
	for n in nodes_inds:
		cursor.execute("SELECT i,start,end "
					   "FROM sets_for_nodes "
					   "WHERE nodes_ind={n} ORDER BY i,start".format(n=n))
		cur_set = []
		new_set = IntervalSet()
		for r in cursor.fetchall():
			intvl = Interval(r['start'],r['end'])
			cur_set.append( (int(r['i']),intvl) )
			new_set.add(intvl)
		if len(new_set) < len(cur_set):
			new_q_sets = defaultdict(lambda: [])
			last_q_index = 0
			for db_i,r in cur_set:
				for q_index,q in enumerate(new_set[last_q_index:]):
					if r.overlaps(q):
						new_q_sets[q_index+last_q_index].append(db_i)
						last_q_index = q_index
						break
			
			# ok, so now, for each q in new_q_sets, there is at least 1 db_i
			# so we use the 1st-entry db_i to update its start,end
			# and set all the other entries to the first db_i, then delete them
			for q_index,db_i_s in new_q_sets.iteritems():
				if len(db_i_s) == 1: continue
				q = new_set[q_index]
				picked_db_i = db_i_s[0]
				others = db_i_s[1:]

				update_node_args.append( (q.lower_bound, q.upper_bound, picked_db_i) )
				update_parsed_args += [(picked_db_i,x) for x in others]
				del_node_args += others

			if len(update_node_args) > 100:
				print >> sys.stderr, "updating......"
				start_t = time.time()
				cursor.executemany("UPDATE sets_for_nodes SET start=%s,end=%s WHERE i=%s", \
						update_node_args)
				cursor.executemany("UPDATE IGNORE parsed SET i1=%s WHERE i1=%s", \
					update_parsed_args)
				cursor.executemany("UPDATE IGNORE parsed SET i2=%s WHERE i2=%s", \
					update_parsed_args)
				cursor.executemany("DELETE FROM parsed WHERE i1=%s", \
					del_node_args)
				cursor.executemany("DELETE FROM parsed WHERE i2=%s", \
					del_node_args)
				cursor.executemany("DELETE FROM sets_for_nodes WHERE i=%s", \
					del_node_args)
				
				update_node_args = []
				update_parsed_args = []
				del_node_args = []
				print >> sys.stderr, "time for 100: ",(time.time()-start_t)
		else:
			print >> sys.stderr, "no update necessary for node index {0}".format(n)

	if len(update_node_args) > 0:
		print >> sys.stderr, "running the remainder sql cmds...."
		cursor.executemany("update sets_for_nodes set start=%s,end=%s WHERE i=%s", \
				update_node_args)
		cursor.executemany("UPDATE IGNORE parsed set i1=%s WHERE i1=%s", \
				update_parsed_args)
		cursor.executemany("UPDATE IGNORE parsed SET i2=%s WHERE i2=%s", \
				update_parsed_args)
		cursor.executemany("DELETE FROM parsed WHERE i1=%s", \
				del_node_args)
		cursor.executemany("DELETE FROM parsed WHERE i2=%s", \
				del_node_args)
		cursor.executemany("DELETE FROM sets_for_nodes WHERE i=%s", \
				del_node_args)

	cursor.close()
	conn.close()

def export_to_db(G, nodes_to_index, i, output_prefix):
	"""
	Given the homology graph G, nodes_to_index, starting i ('i' field in sets_for_nodes)
	Outputs to <output_prefix>.sets_for_nodes and <output_prefix>.parsed
	"""
	sets_for_nodes = {}
	for id,n in nodes_to_index.d.iteritems():
		# since nodes_to_index contains all IGR heads, it's possible the
		# blast files we just read didn't contain some of them
		if not G.has_node(n):
				continue
		print >> sys.stderr, n
		i_set = IntervalSet()
		for e in G.edges_iter(n, data=True):
			# e = (n, other_node_index, BlastEdge instance)
			i_set.update(e[2].lines_of(n))
		sets_for_nodes[n] = i_set

	node_sets_i_map = {}
	with open(output_prefix+'.sets_for_nodes', 'w') as f:
		for n,intervals in sets_for_nodes.iteritems():
			for count,r in enumerate(intervals):
				f.write("{0}\t{1}\t{2}\t{3}\n".format(i,n,r.lower_bound,r.upper_bound))
				node_sets_i_map[(n,count)] = i
				i += 1

	with open(output_prefix+'.parsed', 'w') as f:
		for e in G.edges_iter(data=True):
			node_ind1 = min(e[0],e[1])
			node_ind2 = max(e[0],e[1])
			i1s = set()
			i2s = set()
			for line1 in e[2].lines_of(node_ind1):
				# find the sets_to_nodes index for i1's interval
				for count,r in enumerate(sets_for_nodes[node_ind1]):
					if r.overlaps(line1):
						i1s.add(node_sets_i_map[(node_ind1, count)])
						break
			for line2 in e[2].lines_of(node_ind2):
				# find the sets_to_nodes index for i2's interval
				for count,r in enumerate(sets_for_nodes[node_ind2]):
					if r.overlaps(line2):
						i2s.add(node_sets_i_map[(node_ind2, count)])
						break
			opp = e[2].opposite_strand*1		
			for x,y in itertools.product(i1s, i2s):
				f.write("{0}\t{1}\t{2}\n".format(x,y,opp))

	return i

class BlastEdge:
	def __init__(self, i1, i2, intvl1, intvl2, score, opposite_strand):
		self.lines = {i1: IntervalSet([Interval(*intvl1)]), \
			      i2: IntervalSet([Interval(*intvl2)])}
		self.score = score
		self.muddy = False
		self.opposite_strand = opposite_strand
	
	def __str__(self):
		it = self.lines.iteritems()
		k1,v1 = it.next()
		k2,v2 = it.next()
		return "{0}:{1}\t{2}:{3}\t{4}({5})\t{6}".format(k1, v1, k2, v2, \
				self.score, self.muddy, self.opposite_strand)

	def lines_of(self, i):	
		return self.lines[i]
		
	def add(self, i1, i2, r1, r2, score):
		self.lines[i1].add(Interval(*r1))
		self.lines[i2].add(Interval(*r2))
		if score > self.score:
			self.score = score
			self.muddy = True

def step1_process_blast(blast_output, score_cutoff, nodes_to_index, G, G_nodes_limit=2.5*10**6, program):
	"""
	p.s. now using networkx 0.99, which stupidly doesn't ban selfloops!
	"""
	if program == 'WU':
		blast_parse_func = parseWUBLASTline
	else:
		raise Exception, "temporarily does not support non-WU BLAST output!"

	with open(blast_output) as f:
		line_count = 0
		for line in f:
			line_count += 1
			if G.number_of_nodes() >= G_nodes_limit:
				print >> sys.stderr, "EXITING AT LINE COUNT {0} DUE TO EXCEEDING GRAPH LIMIT".format(line_count)
				return

			d = blast_parse_func(line)
			# ignore self-hits & low score hits
			if d['id1'] == d['id2'] or d['sprime'] < score_cutoff:
				continue

			i1 = nodes_to_index[d['id1']]
			i2 = nodes_to_index[d['id2']]

			if d['start1'] < d['end1']: x1 = (d['start1'],d['end1'])
			else:                       x1 = (d['end1'],d['start1'])
			if d['start2'] < d['end2']: x2 = (d['start2'],d['end2'])
			else:                       x2 = (d['end2'],d['start2'])
			opposite_strand = (d['strand1']!=d['strand2'])
			# if the edge exists, we just expand it
			# (default Graph behaviour is overwrite same edges)
			if G.has_edge(i1, i2):
				e = G.get_edge(i1, i2)
				if e.opposite_strand == opposite_strand:
#					print >> sys.stderr, "before add it is....{0}".format(e)
					e.add(i1, i2, x1, x2, d['sprime'])
#					print >> sys.stderr, "and we're adding to it {0}:{1}:{2},{3}:{4}:{5}".format(i1,d['id1'],x1,i2,d['id2'],x2)
#					print >> sys.stderr, "the line was ", line.strip()
#					print >> sys.stderr, "after add becomes...{0}".format(e)
#					raw_input()
				elif d['sprime'] > e.score:
					print >> sys.stderr, "uh-oh, swapped strand for {0},{1} from {2} to {3}".format(i1,i2,e.score,d['sprime'])
#					print >> sys.stderr, "before the swap the edge was {0}".format(e)
#					print >> sys.stderr, "and we are swapping it with {0},{1}".format(x1,x2)
#					print >> sys.stderr, "the line was ", line.strip()
#					raw_input()
					G.add_edge(i1, i2, BlastEdge(i1, i2, x1, x2, d['sprime'], opposite_strand))
			else:
				G.add_edge(i1, i2, BlastEdge(i1, i2, x1, x2, d['sprime'], opposite_strand))

if __name__=="__main__":
	import_igr_heads_into_db(sys.argv[1])
