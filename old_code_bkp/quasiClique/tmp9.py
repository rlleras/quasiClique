from cPickle import *
import cluster_steps2 as c2
from miscMySQL import get_conn_Actino, get_dict_cursor
from collections import defaultdict
from operator import itemgetter
from scipy import histogram
"""
	This tmp script is used to do evaluation on the cliques we generated
"""
clique_filename='output/output_cliques/ALLActino_RefSeq25_m30s0_cut35gamma80.cliques.pickle'
# here we keep track of some stats
ncRNA_ids_seen = defaultdict(lambda: 0)
clique_stats = defaultdict(lambda: {'sizes':[],'precisions':[]})
ncRNA_id_to_family = {}
hitQ_sizes = []

conn = get_conn_Actino()
cursor = get_dict_cursor(conn)
with open(clique_filename) as f:
	QQQ = load(f)
	for Q in QQQ:
		if 75 in Q or 163 in Q: continue
		if len(Q) < 5:
			continue
		tally_by_family = defaultdict(lambda: 0)
		for i in Q:
			# if it's not a hit, will return (None,None)
			(ncRNA_id, ncRNA_family) = c2.check_hit(i, cursor)
			tally_by_family[ncRNA_family] += 1
			ncRNA_ids_seen[ncRNA_id] += 1
			ncRNA_id_to_family[ncRNA_id] = ncRNA_family
		# decide the dominant family of this cluster
		tally_by_family = tally_by_family.items()
		tally_by_family.sort(key=itemgetter(1))
		fam,count = tally_by_family[-1]
		if fam is not None:
			print("{0}\t{1}/{2}".format(fam,count,len(Q)))
			hitQ_sizes.append(len(Q))
		clique_stats[fam]['sizes'].append(len(Q))
		clique_stats[fam]['precisions'].append(count*1./len(Q))

# now print the stats
print('####################### stats #########################')
print('--- unique hits per family ---')
unique_hits = defaultdict(lambda: 0)
for k in ncRNA_ids_seen: unique_hits[ncRNA_id_to_family[k]] += 1
for fam,count in unique_hits.iteritems():
	print("{0}\t{1}".format(fam,count))
# build a histogram for # of times the same ncRNA was in some clique
hist = defaultdict(lambda: 0)
for v in ncRNA_ids_seen.itervalues(): hist[v] += 1
print('--- histogram of same ncRNA included in some clique ---')
print("1\t{0}".format(hist[1]))
print("2\t{0}".format(hist[2]))
print("3-5\t{0}".format(hist[3]+hist[4]+hist[5]))
print("6-10\t{0}".format(sum(map(lambda i:hist[i], xrange(6,11)))))
print("10-20\t{0}".format(sum(map(lambda i:hist[i], xrange(11,21)))))
print("> 20\t{0}".format(sum(map(lambda i:hist[i], xrange(21,9999)))))
print('--- stats for different families ---')
print('ncRNA family\t# of cliques\tavg.clique size\tavg. precision')
for fam,stats in clique_stats.iteritems():
	print("{0}\t{1}\t{2:.2f}\t{3:.2f}".format(fam,len(stats['sizes']),sum(stats['sizes'])*1./len(stats['sizes']),sum(stats['precisions'])*1./len(stats['precisions'])))
print('---- cliques with hits sizes histogram ----')
print(histogram(hitQ_sizes,[5,11,21,51,101,201]))
