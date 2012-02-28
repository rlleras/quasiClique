import os,re,sys
from networkx import *
from random import seed, choice, uniform
from bisect import bisect
from time import time
from scipy import sparse
from psyco.classes import *

seed(time())

def quick_print(QQ, nodes):
	for Q in QQ:
		Q = map(lambda x: nodes[x][:nodes[x].find('_')], Q)
		print Q.count('moco'),len(Q)

def construct(S, H, alpha, starting_node):
	Q = [starting_node]
	C = S[starting_node, :].nonzero()[1]
	len_C = len(C)
	while len_C > 0:
#		degs_of_C = map(lambda x: S[x, C].nnz, C)
		degs_of_C = H[:, C].sum(axis=1).getA1() # this is much after than the above line!
		min_deg_C = min(degs_of_C)
		max_deg_C = max(degs_of_C)
		RCL_threshold = min_deg_C + alpha*(max_deg_C - min_deg_C)
		# we've randomly picked u from RCL
		RCL = filter(lambda i: degs_of_C[i] > RCL_threshold, xrange(len_C))
		# means there are now fitting RCLs!
		if len(RCL) == 0:
			break
		u = C[choice(RCL)]
		print >> sys.stderr, "picking ",u
		Q.append(u)
		C = C[S[u, C].nonzero()[1]]
		len_C = len(C)
	return Q

def local(H, Q, gamma):
	n = H.shape[0]
	len_Q = len(Q)
	h = H[:,Q]
	h_summed = h.sum(axis=1)
	h_summed2 = h_summed.getA1() # using this speeds on the filter function a LOT....@_@
	gamma_threshold = gamma*(len_Q-1)
	cand = filter(lambda i: h_summed2[i]>=gamma_threshold and i not in Q, xrange(n))
	len_cand = len(cand)
	if len_cand == 0:
		return False
	x = h[cand]
	y = x * x.transpose()
	y.setdiag([0]*len_cand)
	y = y.toarray()
	Q_index_set = set(range(len_Q))
	for v in xrange(len_cand):
		u = y[v,:].argmax()
		if y[v,u] >= gamma*(len_Q-1):
			# this is a good (2,1)-exchange pair!
			# note that w is the "index" in Q, Q[w] is the real thing we're removing
			# similarly, we're adding in cand[u] and cand[v]
			w = Q_index_set.difference(set(x[v,:].nonzero()[1]))
			if len(w) > 0:
				# Q = Q U {u,v}\{w}
				Q.pop(w.pop())
			Q += [cand[u], cand[v]]
			return True
	return False	

def grasp(S, H, gamma, maxitr, given_starting_node=None):
	N = H.shape[0]
	H_deg = H.sum(axis=1).getA1().tolist()
	start_t = time()
	x = filter(lambda i: H_deg[i]>=5, xrange(N))
	print >> sys.stderr, "handling H_deg candidates took {0} secs".format(time()-start_t)
	if len(x) == 0: return []

	bestQ = []
	# pick a starting node unless given
	if given_starting_node is None:
#		starting_node = choice(x)
#		x.remove(starting_node)
		x = zip(H_deg, range(N))
		x.sort()
		starting_node = x[-1][1]
		x = x[:-1]
	else:
		starting_node = given_starting_node
#	curitr = min(H_deg[starting_node], maxitr)
#	print >> sys.stderr, "running with starting node {0} for {1} iterations....".format(starting_node, curitr)
	for k in xrange(maxitr):
		# randomly pick alpha uniformly from [0.1,0.9]
		alpha = uniform(0.1, 0.9)
		print >> sys.stderr, "picked starting node {0} with alpha {1}".format(starting_node, alpha)
		Q = construct(S, H, alpha, starting_node)
		if len(Q) <= 1:
			# no valid local exchange can be done...just give up this round
			if given_starting_node is None and len(x) > 0:
				starting_node = x[-1][1]
				x = x[:-1]
			else:
				return []
			print >> sys.stderr, "CHANGING STARTING NODE"
			continue
		print >> sys.stderr, "before local exchange, size is {0}".format(len(Q))
		while local(H, Q, gamma): pass
		print >> sys.stderr, "max clique with {0} as starting node has size {1}".format(starting_node, len(Q))
		if len(Q) > len(bestQ): bestQ = list(Q)
	return bestQ		

def remove_low_deg_rows(H, valid_rows):
	i = 0
	while i < len(valid_rows):
		if H[valid_rows[i],:].nnz <= 3:
			valid_rows.pop(i)
		else:
			i += 1

def calc_clique_fuzziness(G, Q):
	result = 0
	for i,x in enumerate(Q):
		for y in Q[(i+1):]:
			if G.has_edge(x,y):
				result += 1
	return result
def get_cliques(S, H, QQ, gamma, maxitr, valid_rows):
	validS = S[valid_rows,:]
	validS = validS[:,valid_rows]
	validH = H[valid_rows,:]
	validH = validH[:,valid_rows]
	Q = grasp(validS,validH,gamma,maxitr)
#	H = H.tolil()
#	S = S.tolil()
#	for i in Q:
#		H[i, Q] = 0
#		S[i, Q] = 0
#		# we can also just completely zero the rows that have nnz lower than 5
#		# (the threshold 5 is subject to change)
#		if H[i, :].nnz <= 5: H[i, :] = 0
#		if S[i, :].nnz <= 5: S[i, :] = 0
#	H = H.tocsc()
#	S = S.tocsr()

	# attempt 1: use "valid_rows" to maintain what's left
	Q = map(lambda x: valid_rows[x], Q)
	for i in Q:
		valid_rows.remove(i)
	if len(Q) >= 5:
		QQ.append(Q)
		return True
	else:
		return False

def convert_graph_connectivity_to_sparse(G, nodes):
	"""
		Given a networkx graph, return sparse adjacency matrix S and H
		S and H are different in that S's entires contain edge weights
		(if there are multiple edges, behavior is overwrite),
		and H just has a 1 for every non-zero entry.

		The edge data right now is ((strand1, start1, end1),(strand2, start2, end2), score)
	"""
	n = len(nodes)
	S = sparse.lil_matrix((n,n))
	H = sparse.lil_matrix((n,n))
	nodes_to_index = dict(zip(nodes,range(n)))
	for e in G.edges_iter(data=True):
		i = nodes_to_index[e[0]]
		j = nodes_to_index[e[1]]
		try:
			w = e[2][2]
		except:
			w = e[2]
		S[i,j] = w
		S[j,i] = w
		H[i,j] = 1
		H[j,i] = 1
	# we do a lot of column-slicing, so convert to CSC for efficiency	
	S = S.tocsc()	
	H = H.tocsr()

	return S,H

def convert_to_sparse(Neighbors):
	n = len(Neighbors)
	s = sparse.lil_matrix((n,n))
	for i in xrange(n):
		for j in xrange(n):
			if Neighbors[i][j] == '1':
				s[i,j] = 1
	return s

class CliqueSubList:
	Neighbors = None
	IDmap = None

	def __init__(self, **kwargs):
		"""
			[input arguments]
			-- shared: a list of IDs that's the k-1 shared vertices
			-- kth   : a list of IDs that's the k-th vertex for each clique in this sublist
			-- common_neighbors: a BITSTRING representation of common neighbors, ex: '1101'
					     where a '1' at the i-th index indicates that i is connected
					     to all vertices in shared

			pre: global variable Neighbors and IDmap must already be filled out, so the IDs
			     (usually ints) can be mapped back through IDmap, and the Neighbors (also 
			     bit strings) can be retrieved using the IDs.
		"""
		self.shared = kwargs['shared']
		self.kth = kwargs['kth']
		self.common_neighbors = kwargs['common_neighbors']
		self.k = len(self.shared) + 1

	def __str__(self):
		return "shared: " + " ".join(map(lambda k: str(CliqueSubList.IDmap[k]),self.shared)) + "\n" + \
			str(self.k) + "-th: " + " ".join(map(lambda k: str(CliqueSubList.IDmap[k]),self.kth)) + "\n" + \
		       "common: " + self.common_neighbors

	def str_maximal_clique(self, v, u):
		"""
			Print the maximal clique (self.shared, v, u)
		"""
		return str(len(self.shared)+2)," ".join(map(lambda k: str(CliqueSubList.IDmap[k]), self.shared+[v,u]))

	@staticmethod
	def setup_from_Graph(G):
		CliqueSubList.IDmap = {}
		CliqueSubList.Neighbors = {}
		CliqueList = []

		nodes = G.nodes()
		nodes.sort() # must do!
		len_nodes = len(nodes)
#		return nodes,CliqueSubList
		for i,n in enumerate(nodes):
			print >> sys.stderr, "setting up node {0}/{1}...".format(i,len_nodes)
			CliqueSubList.IDmap[i] = n
			bitarr = ''
			kth = []
			# calculate all reachable notes with path length <= 2
			reachables = path.single_source_shortest_path_length(G, n, 2)
			for j in xrange(len_nodes):
				# THIS IS THE NON-FUZZY VERSION
				if G.has_edge(n,nodes[j]):
					bitarr += '1'
					if j > i:
						kth.append(j)
				else: bitarr += '0'
				# THIS IS THE FUZZY VERSION
#				if nodes[j] in reachables:
#					if reachables[nodes[j]] == 1:
#						# directly connected
#						bitarr += '1'
#						if j > i:
#							kth.append(j)
#					elif reachables[nodes[j]] == 2:
#						# has a path length of 2
#						tmp = filter(lambda x: G.has_edge(x,nodes[j]), G.neighbors_iter(n))
#						if len(tmp) >= 10:
#							print >> sys.stderr, "adding fuzzy k-th {0}...".format(nodes[j])
#							if j > i:
#								kth.append(j)
#							bitarr += '1'
#						else:
#							bitarr += '0'
#					else:
#						bitarr += '0'
#				else:
#					bitarr += '0'
					
			CliqueSubList.Neighbors[i] = bitarr
			if len(kth) > 0:
				CliqueList.append(CliqueSubList(shared=[i], kth=kth, common_neighbors=bitarr))
		return CliqueList

	def find_cliques(self):
		newList = []
		maximals = []
		for i,v in enumerate(self.kth[:-1]):
			new_common = int(self.common_neighbors,2) & int(CliqueSubList.Neighbors[v],2)
			new_kth = []
			for u in self.kth[(i+1):]:
				if CliqueSubList.Neighbors[u][v]=='1': # exists edge (v,u)
					if new_common & int(CliqueSubList.Neighbors[u],2) != 0:
						# (shared,v,u) not maximal since they still share at least 1 vertex
						# that's not in {shared,v,u}
						new_kth.append(u)
					else: # (shared,v,u) is a maximal clique
						maximals.append((v,u))
			if len(new_kth) > 0:
				newList.append(CliqueSubList(shared=self.shared+[v], kth=new_kth, common_neighbors=bin(new_common)[2:]))
		return (newList, maximals)

def maximal_cliques(G):
	"""
		Returns (results, IDmap)
		  where IDmap is dict, maps ID (probably int) --> real node name from graph G
		        results is dict, key is the size of the maximal clique,
                        and value is a list of 2-tuples (shared,[(v,u),(v',u'),(v'',u'')....])
			which means that (shared,v,u), (shared,v',u'), (shared,v'',u'')... are 
		     	maximal cliques
	"""
	from collections import defaultdict
	results = defaultdict(lambda: [])
	tmp = CliqueSubList.setup_from_Graph(G)
	new_tmp = []
	while 1:
		for t in tmp:
			(new_t, maximals) = t.find_cliques()
			new_tmp += new_t
			if len(maximals) > 0:
				results[t.k+1].append((t.shared,maximals))
		if len(new_tmp) == 0: break
		(tmp, new_tmp) = (new_tmp, [])
	return (results,CliqueSubList.IDmap)

def test():
	G = Graph()
	G.add_edge(('b','c'))
	G.add_edge(('b','d'))
	G.add_edge(('b','e'))
	G.add_edge(('b','f'))
	G.add_edge(('b','a'))
	G.add_edge(('a','c'))
	G.add_edge(('a','d'))
	G.add_edge(('a','e'))
	G.add_edge(('c','d'))
	G.add_edge(('c','e'))
	G.add_edge(('c','f'))
	G.add_edge(('c','g'))
	G.add_edge(('d','e'))
	G.add_edge(('d','g'))
	G.add_edge(('e','g'))
	G.add_edge(('f','g'))
	print maximal_cliques(G)
	sys.exit(-1)

def run_test(G):	
	tmp = CliqueSubList.setup_from_Graph(G)
	new_tmp = []
	total_maximals = []
	while 1:
		for t in tmp:
			#print '----------'
			#print t
			(new_t, maximals) = t.find_cliques()
			if len(maximals) > 0:
				total_maximals += [(t.shared,maximals)]
#			for v,u in maximals:
#				print t.str_maximal_clique(v,u)
			new_tmp += new_t
		if len(new_tmp) == 0: break
		(tmp, new_tmp) = (new_tmp, [])
#	return tmp
	return t,total_maximals	

#def fuzzy_combine_maximals(T,total_maximals):
#	for shared,v_u_s in total_maximals:
		

if __name__=="__main__":
	test()
