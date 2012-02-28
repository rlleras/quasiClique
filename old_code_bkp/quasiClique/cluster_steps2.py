from cluster_steps_settings import *
from random import choice, random
import pClique as p

"""
=======================================================
Step 2 of Quasi-Clique Algorithm
=======================================================
Runs the Quais-Clique Algorithm
"""
def pre_clique_setup(conn, cursor, min_deg):
	"""
	(1) Remove low-degree nodes from table
	(2) Create seeds table
	"""
	print >> sys.stderr, "calculate node degrees...."
	degs_of_nodes = calc_degs_from_db(conn)
	more = remove_low_degs_from_db(degs_of_nodes, min_deg)
	while more:
		print >> sys.stderr, "remove deg < {0} nodes....".format(min_deg)
		more = remove_low_degs_from_db(degs_of_nodes, min_deg)

	print >> sys.stderr, "cleaning and populating seed table...."
	cursor.execute("delete from seeds")
	cursor.execute("insert ignore into seeds select distinct i1 from parsed")
	cursor.execute("insert ignore into seeds select distinct i2 from parsed")

def run_clique_finding(iterations, QQQ, output_picklename, cursor, seed_deg_range=(100,500)):
	"""
	Run clique-finding algorithm for <iterations>, in each iteration
	picking randomly a node from table seeds.
	
	Return True if executed all iterations with no problem,
	False if has trouble fetching a seed (which probably means seed table empty)
	"""
	for blahblah in xrange(iterations):
		cursor.execute("SELECT i FROM seeds ORDER BY rand() LIMIT 1")
		try:
			seed_i = cursor.fetchone()[0]
		except:
			print >> sys.stderr, "COULD NOT FETCH A SEED. STOP"
			return False
		
		if seed_deg_range is not None:
			cursor.execute("SELECT count(*) FROM parsed WHERE i1={0} OR i2={0}".format(seed_i))
			seed_deg = cursor.fetchone()[0]
			if seed_deg < seed_deg_range[0] or seed_deg > seed_deg_range[1]:
				print >> sys.stderr, "seed {0} has deg {1}, pick another one. Rewind iteration count.".format(seed_i, seed_deg)
				blahblah = blahblah - 1
				continue
#			raw_input("starting node {0} has degree {1}. press to continue....".format(seed_i, seed_deg))

		# immediately delete this from seeds
		cursor.execute("DELETE FROM seeds WHERE i={0}".format(seed_i))
		print >> sys.stderr, "growing a graph using seed {0}".format(seed_i)
		(Q0,Q) = clique_with_peel(seed_i, cursor)
		if len(Q) == 0:
			print >> sys.stderr, "the graph sucked. Try another seed"
			continue

		# simply delete all existence of this node from parsed
		# TODO: is this the right thing to do?
		print >> sys.stderr, "deleting nodes of the clique from parsed....."
		#cursor.executemany("DELETE IGNORE FROM parsed where i1=%s and i2=%s", itertools.combinations(Q,2))
		tmp = ",".join( map(str, Q) )
		cursor.execute("DELETE IGNORE from parsed where i1 in ({0}) or i2 in ({0})".format(tmp))
		# deleting all nodes from graph also means they're now useless in seeds, so delete 'em too
		cursor.execute("DELETE IGNORE from seeds where i in ({0})".format(tmp))

		if len(Q) >= CLIQUE_MIN_SIZE:
			QQQ.append((Q0,Q)) #QQQ.append(Q)

		with open(output_picklename, 'w') as handle:
			dump(QQQ, handle)
		print >> sys.stderr, "PICKLED to {0}!!!".format(output_picklename)

	return True

def clique_with_peel(seed_i, cursor):
	THRESHOLD = CLIQUE_MIN_SIZE * QUASI_GAMMA - 1
	G = seed_graph(seed_i, cursor)
	# first time we're lenient
	G.delete_nodes_from( filter(lambda n: G.degree(n)<=THRESHOLD, G.nodes_iter()) )
	print >> sys.stderr, "initial seed graph gets {0} nodes, {1} edges".format(G.number_of_nodes(),G.number_of_edges())

	G_nodes = G.nodes()
	S,H = p.convert_graph_connectivity_to_sparse(G, G_nodes)
	tQ = p.grasp(S, H, QUASI_GAMMA, QUASI_MAXITR, G_nodes.index(seed_i))
	Q = map(lambda x: G_nodes[x], tQ)

	print >> sys.stderr, "before peel clique size {0}, graph {1} nodes, {2} edges".format(len(Q),\
			G.number_of_nodes(), G.number_of_edges())

	if seed_i not in Q:
		print >> sys.stderr, "seed node got kicked out PART 1. NOOOOO!!!"
		return ([],[])  #cHANGE BACK LATER

	peel(G, seed_i, Q, QUASI_GAMMA, cursor)
	print >> sys.stderr, "after peel seed graph gets {0} nodes, {1} edges".format(G.number_of_nodes(),G.number_of_edges())

#	return Q, G # TODO: delete later

	G_nodes = G.nodes()
	S,H = p.convert_graph_connectivity_to_sparse(G, G_nodes)
	print >> sys.stderr, "index of seed i", G_nodes.index(seed_i)
	tQ = p.grasp(S, H, QUASI_GAMMA, QUASI_MAXITR, G_nodes.index(seed_i))
	# TESTING (TODO: DELETE LATER) for quasi=0.8 followed by quasi=0.6
	print >> sys.stderr, "after peel clique size {0}".format(len(tQ))
	Q0 = map(lambda x: G_nodes[x], tQ)
	p.local_extra(H, tQ, 0.6)
	print >> sys.stderr, "after local quasi=0.6 clique size {0}".format(len(tQ))
	
	Q = map(lambda x: G_nodes[x], tQ)

	if seed_i not in Q:
		print >> sys.stderr, "seed node got kicked out PART 2. NOOOOO!!!"
		return ([],[])

	# sanity check, delete later (TODO)
	if len(Q) <= CLIQUE_MIN_SIZE:
		return (Q0,Q) # this line is just for avoiding error at the sanity check, delete later
	threshold = 0.6 * len(Q) # CHANGE THIS BACK LATER #threshold = QUASI_GAMMA * len(Q)
	for q in Q:
		assert sum( map(lambda x: G.has_edge(q,x), Q) ) >= threshold

	return (Q0, Q) # CHANGE BACK TO return Q later

def peel(G, seed_i, Q, gamma, cursor):
	# now that we now the smallest possible clique has size |Q|
	# add in length-2 neighbors with degree >= gamma*(|Q|-1)
	len_Q = len(Q)
	threshold = gamma*(len_Q-1)
	
	tmp = ",".join( map(str, G.neighbors_iter(seed_i)) )
	cursor.execute("SELECT i1,i2 FROM parsed WHERE i1 IN ({0}) OR i2 IN ({0})".format(tmp))
	G.add_edges_from( cursor.fetchall() )
	print >> sys.stderr, "after adding direct neighbor edges there are {0} nodes, {1} edges".format(\
			G.number_of_nodes(), G.number_of_edges())

#	to_delete = filter(lambda n: G.has_edge(seed_i, n) and G.degree(n) < threshold, G.nodes_iter())
#	while len(to_delete) > 0:
#		G.delete_nodes_from( to_delete )
#		to_delete = filter(lambda n: G.has_edge(seed_i, n) and G.degree(n) < threshold, G.nodes_iter())
#	print >> sys.stderr, "after removing low-deg direct neighbors there are {0} nodes {1} edges".format(\
#			G.number_of_nodes(), G.number_of_edges())

	return True


def remove_low_degs_from_db(degs_of_node, min_deg):
	"""
	Careful using this, as we're deleting entries from DB!!
	As it deletes nodes, it updates degs_of_node.
	But NOTE that this is ONLY ONE ROUND of deletion.
	You must repeatedly call this to iteratively remove low-deg nodes.

	Returns True if there were stuff to delete, False otherwise.
	"""
	with CONN_FUNC() as cursor:
		to_delete = filter(lambda i: degs_of_node[i]<=min_deg, degs_of_node)
		if len(to_delete)==0:
			return False
		for i in to_delete:
			print >> sys.stderr, "deleting ",i
			cursor.execute("DELETE FROM seeds WHERE i={0}".format(i))
			cursor.execute("SELECT i1,i2 FROM parsed WHERE i1={0} OR i2={0}".format(i))
			for r in cursor.fetchall():
				degs_of_node[r[0]] -= 1
				degs_of_node[r[1]] -= 1
			cursor.execute("DELETE FROM parsed WHERE i1={0} OR i2={0}".format(i))
			del degs_of_node[i]
	return True			
	

def check_hit(i, cursor=None):
	"""
		Given <i>, piece up its accession #, start, end,
		and call get_ribo1 which will return (<ncRNA_id>,<ncRNA_family>)
		if it is a hit or (None,None) if not a hit
	"""
	from miscRibo import get_ribo1
	del_it = False
	if cursor is None:
		conn = CONN_FUNC()
		cursor = get_dict_cursor(conn)
		del_it = True
	cursor.execute("SELECT n.id,s.start,s.end \
					FROM sets_for_nodes AS s \
					LEFT JOIN nodes_to_index AS n \
					ON (s.nodes_ind=n.ind) WHERE i={0}".format(i))
	r = cursor.fetchone()
	(acc,junk),strand,start,end = parsed_accID(r['id'],True,r['start'],r['end'])
	if del_it:
		conn.close()
	return get_ribo1(acc,start,end)


#def db_remove_seeds_not_a_hit(prob_to_keep, low, high, cursor):
#	"""
#		For entries in table seeds that are not a ncRNA hit,
#		only keep them with a probability of <prob_to_keep>
#	"""
#	cursor.execute("SELECT i FROM seeds WHERE i>={0} AND i<={1}".format(low, high))
#	seeds = map(lambda x: x['i'],cursor.fetchall())
#	for i in seeds:
#		if check_hit(i, cursor)[0] is None:
#			if random() > prob_to_keep:
#				cursor.execute("DELETE FROM seeds WHERE i={0}".format(i))
#		else:
#			print >> sys.stderr, "is a hit:", i
#
#
#def convert_i_to_ids(Q):
#	"""
#		Given a list of i in sets_of_nodes,
#		map them back to a list of ('acc/start-end',loc_start,loc_end)s...
#	"""
#	with CONN_FUNC as conn: 
#		cursor = conn.cursor()
#		result = []
#		for q in Q:
#			cursor.execute("SELECT nodes_ind,start,end FROM sets_for_nodes WHERE i={0}".format(q))
#			nodes_ind,start,end = cursor.fetchone()
#			cursor.execute("SELECT id FROM nodes_to_index WHERE ind={0}".format(nodes_ind))
#			result.append((cursor.fetchone()[0],start,end))
#	return result	
#
#def neighbors_of(i, cursor):
#	friends = []
#	cursor.execute("SELECT i2 from parsed WHERE i1={0}".format(i))
#	friends += map(lambda x: x[0], cursor.fetchall())
#	cursor.execute("SELECT i1 from parsed WHERE i2={0}".format(i))
#	friends += map(lambda x: x[0], cursor.fetchall())
#	return friends

def seed_graph(seed_i, cursor):
	G = Graph()
	friends = ''
	cursor.execute("SELECT i2 FROM parsed WHERE i1={0}".format(seed_i))
	for r in cursor.fetchall():
		G.add_edge(seed_i, r[0])
		friends += str(r[0])+','
	cursor.execute("SELECT i1 FROM parsed WHERE i2={0}".format(seed_i))
	for r in cursor.fetchall():
		G.add_edge(seed_i, r[0])
		friends += str(r[0])+','
	if len(friends) > 0:	
		cursor.execute("SELECT i1,i2 FROM parsed WHERE i1 in ({0}) and i2 in ({0})".format(friends[:-1]))
		for r in cursor.fetchall():
			G.add_edge(r[0],r[1])
	return G	

def calc_degs_from_db(db):
	degs_of_node = defaultdict(lambda: 0)

	db.query("SELECT i1,i2 FROM parsed")
	r = db.use_result()
	while 1:
		tmp = r.fetch_row(1)
		if len(tmp) == 0:
			break
		degs_of_node[int(tmp[0][0])] += 1
		degs_of_node[int(tmp[0][1])] += 1

	return degs_of_node	

#def get_components(dbname='ALLActino'):
#	"""
#		This should be called after running cluster_steps on WU-BLAST output,
#		then having exported them to DB and removed all redundancy.
#	"""
#	import _mysql
#	from unionfind import UnionFind
#	db = _mysql.connect(host="xmen", port=7777, user="root", passwd="fidelcastro", db="ALLActino")
#
#	db.query("SELECT i from sets_for_nodes")
#	r = db.use_result()
#	uf = UnionFind()
#	inds = []
#	while 1:
#		x = r.fetch_row()
#		if len(x)==0: break
#		inds.append(long(x[0][0]))
#	uf.insert_objects(inds)	
#
#	raw_input('done with inserting....')
#	db.query("SELECT i1,i2 from parsed")
#	r = db.use_result()
#	count = 0
#	while 1:
#		xs = r.fetch_row(maxrows=10000)
#		if len(xs)==0: break
#		#if count > 100000: break
#		count += 10000
#		print >> sys.stderr, "inserting....", count
#		for x in xs:
#			uf.union(long(x[0]), long(x[1]))
#	db.close()
#
#	components = defaultdict(lambda: set())
#	for i in inds:
#		components[uf.find(i)].add(i)
#	return inds,components

#inds,comp = get_components()

if __name__=="__main__":
	from optparse import OptionParser
	import platform
	parser = OptionParser()
	parser.add_option("-p", "--prefix", help="Output pickle will be <prefix>_<machine>.pickle")
	parser.add_option("-n", "--iterations", help="number of iterations to run, [def 5000]", type="int", default=5000)
	parser.add_option("-s", "--seedrange", help="seed degree range [def (100,700)]", default=(100,700))
	(options, args) = parser.parse_args()
	
	conn = CONN_FUNC()
	cursor = conn.cursor()
	QQQ = []
	dummy_i = 1
	machine = platform.node()
	machine = machine[:machine.find('.')]
	while os.path.exists(options.prefix+'_'+machine+'_'+str(dummy_i)+'.pickle'):
		dummy_i += 1
	output_picklename = options.prefix+'_'+machine+'_'+str(dummy_i)+'.pickle'	
	if type(options.seedrange) is str:
		seedrange = eval(options.seedrange)
	else:	
		seedrange = options.seedrange

	os.system("touch " + output_picklename)	
	print("CONNECTION: {conn}".format(conn=CONN_FUNC))
	print("QUASIGAMMA: {0}".format(QUASI_GAMMA))
	print("ITERATIONS: {iter}".format(iter=options.iterations))
	print("OUTPUT    : {0}".format(output_picklename))
	print("SEEDRANGE : {0}".format(seedrange))
	raw_input("--- is this OK?")
	run_clique_finding(options.iterations, QQQ, output_picklename, cursor, seedrange)
	conn.close()


