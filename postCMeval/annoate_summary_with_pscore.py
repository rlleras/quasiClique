import os
import sys
import glob
import csv

def annotate_rank_summary_with_pscore(filename, delimiter=','):
	"""
	Given a rank_cmfinder.pl-output summary file X, 
	create a new file X.pscore_added that has the motifs' pscores appended
	"""
	motif_dir = os.path.dirname(filename)
	in_csv = csv.DictReader(open(filename), delimiter=delimiter) # cmfinder rank summaries are comma-separated
	with open(filename + '.pscore_added', 'w') as out_f:
		if in_csv.fieldnames is None:
			print >> sys.stderr, "file {0} is odd. IGNORE now".format(filename)
			return	
		new_fieldnames = in_csv.fieldnames + ['pscore']
		out_csv = csv.DictWriter(out_f, new_fieldnames, delimiter=delimiter)
		# need to write out the field names
		#out_csv.writeheader()# lol this function only in 2.7 and i have 2.6 Orz 
		out_f.write(delimiter.join(new_fieldnames) + '\n')
		for obj in in_csv:
			motif_full_path = os.path.join(motif_dir, obj['motif'])
			pscore = os.popen("grep \"Total pair posterior\" {0}.pscoreout".format(motif_full_path)).read().strip()
			obj['pscore'] = float( pscore[len('Total pair posterior '):] )
			out_csv.writerow(obj)

if __name__ == "__main__":
	for filename in glob.iglob('motifs/*/*.fna.summary'):
		print >> sys.stderr, "annotating pscore to {0}.pscore_added....".format(filename)
		annotate_rank_summary_with_pscore(filename)
