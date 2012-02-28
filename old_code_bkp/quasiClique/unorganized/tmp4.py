import os,re,sys
from networkx import *
from networkx.readwrite import *
from mlabwrap import mlab
from numpy import array

SMALLEST_SUBGRAPH_SIZE = 5

g = read_gpickle('yybPM8N7Q16R2.WUblast_cut30.step2ed.gpickle')
g.remove_all_selfloops()
g.remove_all_multiedges()# this is rough, but leave this this way now
g.ban_selfloops()

#remove loners
g.delete_nodes_from(filter(lambda n: g.degree(n)==0, g.nodes_iter()))

node_mapping = {}
node_mapping_count = 0
f = open('filename.txt','w')
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
	#f.write("{0}\t{1}\t{2}\n".format(i,j,e[2][2]))
	f.write("{0}\t{1}\t1\n".format(i,j))
f.close()	

mlab.readin('filename.txt','miska')

HCS = []
todo = [array(range(1,node_mapping_count))]
bad = []

while len(todo) > 0:
	cur = todo.pop()
	print >> sys.stderr, "running matlab on size {0}.....".format(len(cur))
	answer = mlab.test('miska',cur)
	if answer[0][0]==-1:
		print >> sys.stderr, "bad bunch of size {0}!".format(len(cur))
		bad.append(list(cur))
#	elif answer[0][0] > len(cur)*.9: # is a highly connected subgraph
#		HCS.append(list(cur))
#		print >> sys.stderr, "cut {0} is HCS".format(len(cur))
	else:
		set1 = answer[0][1:]
		set2 = array(list(set(list(cur)).difference(set(list(set1)))))
		if len(set1) > SMALLEST_SUBGRAPH_SIZE:
			todo.append(set1)
		else:
			HCS.append(list(set1))
		if len(set2) > SMALLEST_SUBGRAPH_SIZE:	
			todo.append(set2)
		else:
			HCS.append(list(set2))
		print >> sys.stderr, "{0} split into {1},{2}".format(len(cur),len(set1),len(set2))
	# sanity check
	try:
		assert sum(map(len, HCS)) + sum(map(len, todo)) + sum(map(len, bad)) == node_mapping_count-1
	except:
		print "that was: ",list(cur)
		sys.exit(-1)
#	raw_input()	

reverse_node_mapping = {}
for k,v in node_mapping.iteritems():
	reverse_node_mapping[v] = k

for s in HCS:
	for x in s:
		print reverse_node_mapping[int(x)],
	print ''
	
