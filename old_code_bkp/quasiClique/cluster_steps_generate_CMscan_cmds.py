from cPickle import *
import os,re,sys,fnmatch
from cluster_steps4 import MotifRankInfo

# ---------- parse input options ----------------- #
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-d", "--dir", help="Directory containing directories of motifs")
parser.add_option("-i", "--input", help="Pickle filename, which should contain a list of MotifRankInfo")
parser.add_option("-s", "--scan_filename", \
		help="Filename of the fasta seqs we're going to scan against "
		"default is Zasha_20090819_nonredundant_plus5Xdishuffled.as1chromosome.fna",\
		type="string",
		default='Zasha_20090819_nonredundant_plus5Xdishuffled.as1chromosome.fna')
parser.add_option("-b", "--usebest", help="use the best motif (or else will use every motif)", \
		action="store_true", default=False, dest="usebest")

(options, args) = parser.parse_args()

if options.usebest:
	print >> sys.stderr, "Using best motif"
else:
	print >> sys.stderr, "Using EVERY motif"

with open(options.input) as handle:
	MRIs = load(handle)
	# remove from the list those that have motif family being 'None'
	MRIs = filter(lambda obj: obj.fam is not None, MRIs)

motif_rex = re.compile( '\S+.fna\S+motif\S+[\.h\d+_\d+]+$' )

for obj in MRIs:
	# some hardcoded crap here....*sigh*, right now the motif_filename is stored like:
	# grid01_02_124.fna.1.motif.h1_2.dup_rmed.html (has .dup_rmed. it's used with rankpl)
	motif_dir = obj.motif_filename[:obj.motif_filename.find('.')]
	
	dd = os.path.join( options.dir, motif_dir )
	# right now, look for the NONE dup_rmed version!
	for file in fnmatch.filter(os.listdir(dd), '*motif*'):
#		if file.endswith('.dup_rmed'):
#			continue

		full_filename = os.path.join( dd, file )
		if not open(full_filename).readline().strip().startswith('# STOCKHOLM 1.0'):
			continue

		# TODO: is this the right thing to do?
		best_motif = obj.motif_filename[:obj.motif_filename.rfind('.')]
		if (options.usebest and file == best_motif) or not options.usebest:
			# to use cmcalibrate it would be cmcalibrate {0}.cmbuilded;
			print("cmbuild {0}.cmbuilded {0}; cmcalibrate {0}.cmbuilded;"\
					"cmsearch --fil-T-hmm 10 --tabfile {0}.hascalib_hmmT10.cmsearched_hits {0}.cmbuilded {1}".format(\
							full_filename, options.scan_filename))
