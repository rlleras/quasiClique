from cluster_steps_settings import *
from random import choice, random, uniform
import pClique as p

"""
=======================================================
Step 2 of Quasi-Clique Algorithm
=======================================================
Runs the Quais-Clique Algorithm
"""
def pre_clique_finding_setup():
	dd = calc_degs_from_db()
	while remove_low_degs_from_db(dd):
		pass
	print >> sys.stderr, "clean up seeds table in case there are stuff...."	
	cursor.execute("delete from seeds")
	cursor.exeucte("insert into seeds select distinct i1 from parsed")
	cursor.execute("insert ignore into seeds select distinct i2 from parsed")

def run_clique_finding(iterations, QQQ, output_picklename, cursor, seed_deg_range=(100,700)):
	"""
	Run clique-finding algorithm for <iterations>, in each iteration
	picking randomly a node from table seeds.
	
	Return True if executed all iterations with no problem,
	False if has trouble fetching a seed (which probably means seed table empty)
	"""
	blahblah = 0
	badblood = 0
	while blahblah < iterations and badblood < iterations:
		cursor.execute("SELECT i FROM seeds ORDER BY rand() LIMIT 1")
		try:
			seed_i = cursor.fetchone()[0]
		except:
			print >> sys.stderr, "COULD NOT FETCH A SEED. STOP"
			return False
		# immediately delete this from seeds (TODO: right thing to do here?)
		cursor.execute("DELETE FROM seeds WHERE i={0}".format(seed_i))

		cursor.execute("SELECT count(*) FROM parsed WHERE i1={0} OR i2={0}".format(seed_i))
		if seed_deg_range is not None:
			seed_deg = cursor.fetchone()[0]
			if seed_deg < seed_deg_range[0] or seed_deg > seed_deg_range[1]:
				print >> sys.stderr, "degree", seed_deg, " too low or high, pick another one"
				badblood += 1
				continue
		blahblah += 1
		Q = clique_perfects(seed_i, cursor)
		if len(Q) == 0:
			print >> sys.stderr, "the graph sucked. Try another seed"
			continue
		# NOTE: we remove ALL EDGES between the nodes, not the NODES in this version!!!
		ss = ",".join( map(str, Q) )
		
		cursor.execute("DELETE FROM parsed WHERE i1 IN ({0}) AND i2 IN ({0})".format(ss))
		QQQ.append(Q)
		with open(output_picklename, 'w') as f:
			dump(QQQ, f)
	return True

def clique_perfects(seed_i, cursor):
	G = seed_graph(seed_i, cursor)
	# first find a perfect clique of size >= PERFECT_CLIQUE_MIN_SIZE
	G.delete_nodes_from( filter(lambda n: G.degree(n)<PERFECT_CLIQUE_MIN_SIZE, G.nodes_iter()) )
	if not G.has_node(seed_i):
		return []
	print >> sys.stderr, "initial seed graph gets {0} nodes, {1} edges".format(\
			G.number_of_nodes(),G.number_of_edges())
	G_nodes = G.nodes()
	S,H = p.convert_graph_connectivity_to_sparse(G, G_nodes)
	remembered_perfects = []
	# for PERFECT_MAXITR iterations, store all found perfect cliques containing seed_i
	for itr in xrange(PERFECT_MAXITR):
		alpha = uniform(.1, .9)
		Q = map(lambda x: G_nodes[x], p.construct(S, H, alpha, G_nodes.index(seed_i)))
		Q.sort() # sort it, so we add only distinct cliques
		if seed_i in Q and len(Q) >= PERFECT_CLIQUE_MIN_SIZE\
				and Q not in remembered_perfects:
			remembered_perfects.append(Q)
	if len(remembered_perfects) == 0:
		return []
	remembered_perfects.sort(reverse=True,key=lambda x:len(x))
	return remembered_perfects[0]

def cmfinder_on_clique(Q, fasta_filename, cursor):
	from miscBio import FastaReader
	import tempfile
	FETCH_SQL = "SELECT n.id,s.start,s.end \
				FROM sets_for_nodes s \
				LEFT JOIN nodes_to_index AS n \
				ON (s.nodes_ind=n.ind) WHERE i={i}"
	F = FastaReader(fasta_filename)
	handle = tempfile.NamedTemporaryFile(suffix='.fna',dir='harwood/MOTIF_DIR/PseudomonasCOG',\
			delete=False)
	for q in Q:
		print >> sys.stderr, "fetching info for ", q
		cursor.execute( FETCH_SQL.format(i=q) )
		id,start,end = cursor.fetchone()
		handle.write(">{id}/{start}-{end}\n{s}\n".format(id=id,start=start,end=end,\
				s=F[id].seq[start:(end+1)]))
	print >> sys.stderr, "file is....",handle.name
	handle.close()

def peel(G, seed_i, len_Q, gamma, cursor):
	# now that we now the smallest possible clique has size (|Q|-1)
	# we first remove all neighbors of degree < gamma*(|Q|-1)
	# then add in length-2 neighbors with degree >= gamma*(|Q|-1)
	threshold = gamma*(len_Q-1)
	if G.has_node(seed_i):
		tmp = ",".join( map(str, G.neighbors_iter(seed_i)) )
	else: # wow, seed node was lame! it got kicked.oh well.
		return False
	if len(tmp) > 0:	
		cursor.execute("SELECT i1,i2 FROM parsed WHERE i1 IN ({0}) OR i2 IN ({0})".format(tmp))
		for r in cursor.fetchall(): G.add_edge(r[0],r[1])
	G.delete_nodes_from( filter(lambda n: G.degree(n) < threshold, G.nodes_iter()) )
	return True


def remove_low_degs_from_db(degs_of_node, min_deg=3):
	"""
	Careful using this, as we're deleting entries from DB!!
	As it deletes nodes, it updates degs_of_node.
	But NOTE that this is ONLY ONE ROUND of deletion.
	You must repeatedly call this to iteratively remove low-deg nodes.

	Returns True if there were stuff to delete, False otherwise.
	"""
	with CONN_FUNC() as cursor:
		to_delete = filter(lambda i: degs_of_node[i]<min_deg, degs_of_node)
		if len(to_delete)==0:
			return False
		for i in to_delete:
			print >> sys.stderr, "deleting ",i
			cursor.execute("DELETE FROM sets_for_nodes WHERE i={0}".format(i))
			cursor.execute("SELECT i1,i2 FROM parsed WHERE i1={0} OR i2={0}".format(i))
			for r in cursor.fetchall():
				degs_of_node[r[0]] -= 1
				degs_of_node[r[1]] -= 1
			cursor.execute("DELETE FROM parsed WHERE i1={0} OR i2={0}".format(i))
			del degs_of_node[i]
	return True			

def seed_graph(seed_i, cursor):
	seed_i = long(seed_i)
	G = Graph()
	friends = set()
	cursor.execute("SELECT i2 FROM parsed WHERE i1={0}".format(seed_i))
	for r in cursor.fetchall():
		G.add_edge(seed_i, r[0])
		friends.add(r[0])
	cursor.execute("SELECT i1 FROM parsed WHERE i2={0}".format(seed_i))
	for r in cursor.fetchall():
		G.add_edge(seed_i, r[0])
		friends.add(r[0])
	if len(friends) > 0:	
		friends = ",".join( map(str, friends) )
		cursor.execute("SELECT i1,i2 FROM parsed \
				WHERE i1 in ({0}) and i2 in ({0})".format(friends))
	for r in cursor.fetchall():
		G.add_edge(r[0],r[1])
	return G	

def calc_degs_from_db(cursor=None):
	degs_of_node = defaultdict(lambda: 0)
	close_it = False
	if cursor is None:
		conn = CONN_FUNC()
		cursor = conn.cursor()
		close_it = True
	cursor.execute("SELECT i1,i2 FROM parsed")
	for r in cursor.fetchall():
		degs_of_node[r[0]] += 1
		degs_of_node[r[1]] += 1
	if close_it:
		conn.close()
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
	print("ITERATIONS: {iter}".format(iter=options.iterations))
	print("OUTPUT    : {0}".format(output_picklename))
	print("SEEDRANGE : {0}".format(seedrange))
	raw_input("--- is this OK?")
	run_clique_finding(options.iterations, QQQ, output_picklename, cursor, seedrange)
	conn.close()
