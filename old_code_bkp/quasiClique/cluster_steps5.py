import os,re,sys
from miscCMF import read_motif, furnish_motif
from myAnchor import findPRs, test2
from numpy import array
from cPickle import *

"""
Used for clustering motifs
"""

def compare_motif(motif_filename1, motif_filename2):
	cons1, rf1 = furnish_motif(*read_motif(motif_filename1))
	cons2, rf2 = furnish_motif(*read_motif(motif_filename2))
	# remove gaps and furnish consensus so it's just ( ) .
	pr1 = findPRs(cons1)[0]
	pr2 = findPRs(cons2)[0]
	#print cons1, rf1, pr1
	#print cons2, rf2, pr2
	return test2(rf1,rf2,cons1,cons2,pr1,pr2,'XV','VX',2000,True,['junk1'],['junk2'])

def cluster_motifs(total_ranks_pickle, output_pickle, motif_dir='tmp/motif_dir/'):
	"""
	Input should be a total_ranks pickle outputted from cluster4
	where it's a sorted list of (rank_index,fam,count,motif_filename,data,count_o,motif_size)s
	"""
	from hcluster import linkage, to_tree
	motifs = []
	distance = []
	total_ranks = load(open(total_ranks_pickle))
	for rank_index,fam,count,motif_filename,data,count_o,motif_size in total_ranks:
		if rank_index < 300:
			continue # TODO:delete or change later?
		motifs.append({'file': motif_filename, 'info': "{0}_{1}".format(fam,count)})
	print >> sys.stderr, "comparing {0} motifs.....".format(len(motifs))
	for i in xrange(len(motifs)):
		print >> sys.stderr, "{0}/{1}".format(i,len(motifs))
		for j in xrange(0,i):
			f1 = motifs[i]['file']
			f2 = motifs[j]['file']
			f1 = motif_dir + f1[:f1.find('.fna')] + '/' + f1
			f2 = motif_dir + f2[:f2.find('.fna')] + '/' + f2
			score = compare_motif(f1, f2)['score']
			# the bigger the score (more similar), the smaller their distance
			# so we invert the score to 10000/score achieve that
			# however be careful that score of -9999999 means we should give it
			# a distance of...say...999999?
			if score > 0:
				distance.append( 10000./score )
			else:
				distance.append( 999999. )
	Z = linkage(distance, method='average')
	# delete below later
	f = open(output_pickle,'w')
	dump({'motifs':motifs,'Z':Z}, f)
	f.close()
	# delete above later
	# now Z is the linkage matrix, convert it to newick with proper naming
	root = to_tree(Z)
	print(hcluster_cnode_to_newick(root, motifs))

def hcluster_cnode_to_newick(root, mapping):
	if root.is_leaf():
		return	mapping[root.get_id()]['info']
	else:
		return '(' + hcluster_cnode_to_newick(root.get_left(), mapping) + ',' + \
				hcluster_cnode_to_newick(root.get_right(), mapping) + ')'
	
if __name__=="__main__":
	cluster_motifs('firm_ALL.total_ranks.pickle','firm_ALL.rank300_motif_sims.pickle')
