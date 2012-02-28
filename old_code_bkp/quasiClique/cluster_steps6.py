from cluster_steps_settings import *
from miscRibo import get_ribo1
from miscParses import parsed_accID, calc_overlap
import miscInfernal
from cluster_steps4 import MotifRankInfo
from math import log
"""
=======================================================
Step 6 of Quasi-Clique Algorithm
=======================================================
Used to evaluate CMscans
"""
def make_eval_func_by_mapping_pickle(mapping_pickle_filename):
	"""
	If the scanning fasta file was just one large concatenated chromosome,
	then we need this mapping pickle which is a list of (0-based start, 0-based end, fam)

	Returns a function that takes input arguments id, (best,loc_start,loc_end)
	and returns rb_fam (which can be 'random' to indicate a scan hit to a none-ncRNA region)
	"""
	with open(mapping_pickle_filename) as handle:
		mapping = load( handle )

	def eval_func(id, (best,loc_start,loc_end)):
		if loc_start > loc_end:
			# swap coords if hit is on minus strand
			(loc_start, loc_end) = (loc_end, loc_start)
		ii = bisect(mapping, (loc_start,loc_end))	
		T = defaultdict(lambda: 0)
		# go <--- way
		for j in xrange(min(ii,len(mapping)-1),0,-1):
			if mapping[j][1] < loc_start:
				break
			T[mapping[j][2]] = max(T[mapping[j][2]], calc_overlap(loc_start, loc_end, mapping[j][0], mapping[j][1]))
		# go ---> way
		if ii < len(mapping)-1:
			for j in xrange(ii+1, len(mapping), +1):
				if mapping[j][0] > loc_end:
					break
				T[mapping[j][2]] = max(T[mapping[j][2]], calc_overlap(loc_start, loc_end, mapping[j][0], mapping[j][1]))
		# take the majority family and return it
		T = T.items()
		T.sort(key=itemgetter(1), reverse=True)
		print >> sys.stderr, "best is", T[0], " out of whole length", (loc_end-loc_start+1)
#		raw_input()
		if T[0][1] >= .5*(loc_end-loc_start+1):
			t = T[0][0]
			return t #return t[:t.index('___')]
		else:	
			return 'random___000'#return 'random'

	return eval_func	

def evaluate_CMscaned(dir, scan_fasta_pickle, total_ranks_pickle, best_e_or_score, cutoff,\
		eval_func=None,\
		suffix='.nocalib_hmmT10.cmsearched_hits'):
	"""
	<dir> should be a directory containing files of format <motif_filename>.cmsearched_hits
	which are outputs from Infernal-1.0's cmsearch use --tabfile option.

	<total_ranks_pickle> can either be from rankpl or pscore. Both will be a list of
	cluster_steps4.MotifRankInfo objects.
	"""
	if eval_func is None:
		eval_func = lambda id,(best,loc_start,loc_end): id[:id.index('___')]

	with open(total_ranks_pickle) as handle:
		total_ranks = load(handle)
	with open(scan_fasta_pickle) as handle:
		total_fam_counts = load(handle)
		# some of the total_fam_counts have keys that are like moco-
		# we'll put them as moco...
		for k in total_fam_counts:
			if k is not None and k.endswith('-'):
				total_fam_counts[k[:-1]] += total_fam_counts[k]
				del total_fam_counts[k]

	best_cmscan_by_family = defaultdict(lambda: {'sens':0., 'prec':0., 'motif':None})
	for_scatterplot = {'real':[], 'random':[]}

	# go through each MotifRankInfo obj, if family is not None, then we did a CM scan for it
	# refer to cluster_steps_generate_CMscan_cmds.py to see why this is the case
	for obj in filter(lambda obj: obj.fam is not None, total_ranks):
		# some hardcoded crap here....*sigh*, right now the motif_filename is stored like:
		# grid01_02_124.fna.1.motif.h1_2.dup_rmed.html (has .dup_rmed. it's used with rankpl)
		filename = obj.motif_filename[:obj.motif_filename.rfind('.html')] + suffix
		if not os.path.exists(os.path.join(dir, filename)):
			print("ERROR: file {0} doesn't exist".format(filename))
			continue

		# returns scanned_id --> (best_e_or_score, start, end)
		print >> sys.stderr, "file is {0}, family should be {1}".format(filename, obj.fam)
		cm_results = miscInfernal.read_tab_cmscanned(os.path.join(dir, filename), best_e_or_score, cutoff)

		if len(cm_results) == 0:
			continue
		tally_by_family = defaultdict(lambda: set())
		for id,(best,loc_start,loc_end) in cm_results.iteritems():
#			print >> sys.stderr, id, loc_start, loc_end
#			(_acc,_junk),_strand,_start,_end = parsed_accID(id,version_split=True,loc_start=loc_start,loc_end=loc_end)
#			print >> sys.stderr, "parsed into", _acc, _start, _end
#			rb_id,rb_fam = get_ribo1(_acc, _start, _end)
#			print >> sys.stderr, "rb is", rb_id, rb_fam
#			raw_input()
#			if rb_fam is not None and rb_fam.endswith('-'):
#				rb_fam = rb_fam[:-1]

#			print >> sys.stderr, "evaling", id, best, loc_start, loc_end
			rb_fam_with_id = eval_func(id, (best,loc_start,loc_end)) # in format fam___id
			rb_fam,rb_id = rb_fam_with_id.split('___')
			tally_by_family[rb_fam].add( rb_id )
			if rb_fam == obj.fam:
				for_scatterplot['real'].append( best)
			else:
				for_scatterplot['random'].append( best )
		fam = obj.fam
		for t in tally_by_family:
			tally_by_family[t] = len(tally_by_family[t])
		if fam not in tally_by_family:
			tally_by_family[fam] = 0 # we need to do this otherwise next two lines with fail
		# calculate sensitivity and precision(PPV)
		sens = tally_by_family[fam]*1./total_fam_counts[fam]
		prec = tally_by_family[fam]*1./len(cm_results)
#		print >> sys.stderr, "motif:", os.path.join(dir,filename), "fam: ", fam, "sens: ", sens, "prec: " , prec
		if sens > best_cmscan_by_family[fam]['sens'] or \
			sens==best_cmscan_by_family[fam]['sens'] and prec > best_cmscan_by_family[fam]['prec']:
			best_cmscan_by_family[fam] = {'sens': sens, 'prec': prec, 'motif': filename}

	fams = best_cmscan_by_family.keys()
	fams.sort(key=str.lower)
	# print the final results in latex table format
	for fam in fams:
		print("{fam} & {count} & {sens:.2f} & {prec:.2f} & {motif}\\\\hline".format(fam=fam, \
				count=total_fam_counts[fam], sens=best_cmscan_by_family[fam]['sens'], \
				prec=best_cmscan_by_family[fam]['prec'],\
				motif=best_cmscan_by_family[fam]['motif']))

	# print for_scatterplot
	print(",".join(map(str, for_scatterplot['real'] ) ) )
	print(",".join(map(str, for_scatterplot['random'] ) ))

if __name__=="__main__":
	import psyco
	psyco.full()

	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-d", "--dir", help="directory containing CM scan output (--tabfile format)")
	parser.add_option("-s", "--scanFastaPickle", help="pickle of fam counts of the scanned fasta,"
			"default ../../output/input_fnas/Zasha_20090819_nonredundant.fna.fam_counts.pickle",\
			default='../../output/input_fnas/Zasha_20090819_nonredundant.fna.fam_counts.pickle')
	parser.add_option("-t", "--totalRanksPickle", help="Total Ranks pickle used to generate scans,"
			"default ALLFirm_sdustM8N7Q16R2W3E2_cut35.quasi80iter20size5_ALL.test1.ranked_rankpl.pickle",\
			default='ALLFirm_sdustM8N7Q16R2W3E2_cut35.quasi80iter20size5_ALL.test1.ranked_rankpl.pickle')
	parser.add_option("-c", "--cutoff", type="float", help="E-value or Score cutoff")
	parser.add_option("-m", "--cutoffMethod", help="use either 'e' (E-value) or 'score' for cutoff")
	parser.add_option("-x", "--suffix", help="suffix for cm hit reports, default .nocalib_hmmT10.cmsearched_hits",\
			default='.nocalib_hmmT10.cmsearched_hits')
	(options,args) = parser.parse_args()

	assert options.cutoffMethod in ('e', 'score')

	# TODO: change back later
	eval_func = make_eval_func_by_mapping_pickle('Zasha_20090819_nonredundant_plus5Xdishuffled.as1chromosome.mapping_pickle')
	evaluate_CMscaned(options.dir, options.scanFastaPickle, options.totalRanksPickle, \
			options.cutoffMethod, options.cutoff, eval_func=eval_func, suffix=options.suffix)
