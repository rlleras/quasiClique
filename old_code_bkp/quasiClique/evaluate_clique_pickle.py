import os,re,sys,operator
from cPickle import *
from collections import defaultdict
import cluster_steps4 as c4
import cluster_steps_settings as cset
"""
Evaluates a directory of motifs that has the input .fnas, then outputs in Latex table format.
This is used for writing the paper....
Given the <output_prefix>, outputs two things:
	--- <output_prefix>.evaled_clique_verbose.tex  # a clique-by-clique spit out
	--- <output_prefix>.evaled_clique_summary.tex  # a concise summary
"""

# ------------------- parsing input options -------------------------- #
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--dir", type="string", help="Dir containing directories of motifs. Use either this or -c.")
parser.add_option("-c", "--clique", type="string", help="Clique pickle. Use either this or -d")
parser.add_option("-p", "--phylum", type="string", help="Phylum, [default] Fimircutes", default="Firmicutes")
parser.add_option("-r", "--recoverable", type="string", \
		help="Recoverable pickle for the input IGR, "
		"default ../../output/input_fnas/ALLFirm_RefSeq25_m30s0.IGRs.recoverable_ncRNAs.pickle",\
		default='../../output/input_fnas/ALLFirm_RefSeq25_m30s0.IGRs.recoverable_ncRNAs.pickle')
parser.add_option("-o", "--output_prefix", help="output will be <output_prefix>.evaled_clique_{summary|verbose}.tex")

(options, args) = parser.parse_args()

if not operator.xor(options.dir is None, options.clique is None):
	raise Exception, "ERROR: Only one and exactly one of -d or -c can be used as input!!"

#MOTIF_DIR = '../../MOTIF_DIR/ALLFirm_sdustM8N7Q16R2_quasi60/'
#PHYLUM = 'Firmicutes'
#IGR_RECOVERABLE_PICKLE = '../../output/input_fnas/ALLFirm_RefSeq25_m30s0.IGRs.recoverable_ncRNAs.pickle'
#CLIQUE_PICKLE = './ALLFirm_sdustM8N7Q16R2W3E2_cut35.quasi60iter20size5_ALL.pickle'
#output_prefix = sys.argv[1]

MOTIF_DIR = options.dir
PHYLUM    = options.phylum
IGR_RECOVERABLE_PICKLE = options.recoverable
CLIQUE_PICKLE = options.clique
output_prefix = options.output_prefix

# read the phylum distribution file
#from miscRibo import read_phylum_distribution
#phylum_dist = read_phylum_distribution()
# fake it for None
#phylum_dist[PHYLUM][None] = 100000.
# read the 'recoverable' number of ncRNAs per class
# (recoverable = currently defined as IGR covering at least 50bp or 50%)
with open(IGR_RECOVERABLE_PICKLE) as handle:
	recoverable = load(handle)
recoverable[None] = 100000.

# here we keep track of some stats, key is family
clique_stats = defaultdict(lambda: {'sizes':[],'precisions':[], 'raw_count_set':set()})
# used for writing to verbose, key is family, val is list of (clique_size, count, clique_name)
spit_out = defaultdict(lambda: [])

if MOTIF_DIR is None:
	with open(options.clique) as handle:
		QQQ = load( handle )
	with cset.CONN_FUNC() as cursor:
		for index,Q in enumerate(QQQ):
			print >> sys.stderr, "doing {0}-th clique.....".format(index)
			fam, count, clique_size, ncRNA_ids_hit = c4.eval_clique(Q[1], cursor)
			if fam is not None:
				spit_out[fam].append( (clique_size, count, Q) )
			clique_stats[fam]['sizes'].append(clique_size)
			clique_stats[fam]['precisions'].append(count*1./clique_size)
			clique_stats[fam]['raw_count_set'] = clique_stats[fam]['raw_count_set'].union( ncRNA_ids_hit )
else:
	for d in os.listdir(MOTIF_DIR):
		dd = os.path.join(MOTIF_DIR, d)
		if os.path.isdir(dd):
			print >> sys.stderr, "doing directory.....", d
			fam, count, clique_size, ncRNA_ids_hit = c4.eval_original_fna( os.path.join(dd, d+'.fna') )
			if fam is not None:
				spit_out[fam].append( (clique_size, count, d) )
			clique_stats[fam]['sizes'].append(clique_size)
			clique_stats[fam]['precisions'].append(count*1./clique_size)
			clique_stats[fam]['raw_count_set'] = clique_stats[fam]['raw_count_set'].union( ncRNA_ids_hit )

fams = spit_out.keys()
fams.sort(key=str.lower) # case-insensitive sort
# now write the spit out
with open(output_prefix + '.evaled_clique_verbose.pickle', 'w') as f:
	dump(dict(spit_out), f)

with open(output_prefix + '.evaled_clique_verbose.tex', 'w') as f:
	f.write("""
		\\begin{table}[!htpb]
		\\begin{flushleft}
		\\caption{\\bf{Table Caption}}
		\\end{flushleft}
		\\begin{tabular*}{0.30\\textwidth}{@{\\extracolsep{\\fill}}|l|r|}
		\\hline
		Class & Specificity \\\\\\hline
		\\hline
		""")
	for fam in fams:
		if fam.find('_') >= 0:
			fam1 = fam.replace('_','$\\_$')
		else:
			fam1 = fam
		spit_out[fam].sort()
		for clique_size,count,d in spit_out[fam]:
			f.write("{0} & {1}/{2}\\\\\\hline\n".format(fam1, count, clique_size))
	f.write("""
			\\hline
			\\end{tabular*}
			\\end{table}
			""")

with open(output_prefix + '.evaled_clique_summary.tex', 'w') as f:
	f.write("""\\begin{table}[!htpb]
		\\begin{flushleft}
		\\caption{\\bf{TITLE HERE}}
		\\end{flushleft}
		\\begin{tabular*}{\\hsize}{@{\\extracolsep{\\fill}}|l|r|r|r|r|r|}
		\\hline
		Class & Recoverable & \\# of cliques & avg. clique size & avg. precision & total sensitivity\\\\
		\\hline
		""")

	fams = fams + [None]
	for fam in fams:
		# we have to massage fam so it doesn't choke latex
		if fam is not None and fam.find('_') >= 0:
			fam1 = fam.replace('_','$\\_$')
		else:
			fam1 = fam
		stats = clique_stats[fam]
		avg_size = int(round(sum(stats['sizes'])*1./len(stats['sizes'])))
		avg_prec = "{0:.2f}".format( round(sum(stats['precisions'])*1./len(stats['precisions']), 2) )
		if fam is not None and fam.endswith('--'):
			r_fam = fam[:-2]
		elif fam is not None and fam.endswith('-'):
			r_fam = fam[:-1]
		else:
			r_fam = fam
		total_sens = "{0:.2f}".format( len(stats['raw_count_set'])*1./recoverable[r_fam] ) 
		f.write("""{family} & {recoverable_count} &
		{num_cliques} & 
		{avg_size} & 
		{avg_prec} & 
		{total_sens}\\\\\\hline
		""".format(family=fam1, \
			recoverable_count=recoverable[r_fam],\
			num_cliques=len(stats['sizes']), avg_size=avg_size,\
			avg_prec=avg_prec, total_sens=total_sens))

	f.write("""\\hline
		\\end{tabular*}
		\\end{table}
		""")
