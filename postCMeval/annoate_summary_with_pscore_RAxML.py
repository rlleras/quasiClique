import os
import sys
import glob
import csv
"""
2010/12/06 
 Now adds two additional fields: pscore_RAxML and fixTreeEstBranch.pscoreout_RAxML
"""

def annotate_rank_summary_with_pscore_RAxML(filename, delimiter=','):
	"""
	Given a .pscore_added file which is already augmented with normal pscore (using pfold's tree)
	Create another file called .pscore_added.RAxML_added with has an additional field of "pscore_RAxML"
	"""
	motif_dir = os.path.dirname(filename)
	in_csv = csv.DictReader(open(filename), delimiter=delimiter) # cmfinder rank summaries are comma-separated
	with open(filename + '.RAxML_added', 'w') as out_f:
		if in_csv.fieldnames is None:
			print >> sys.stderr, "file {0} is odd. IGNORE now".format(filename)
			return	
		pscore_fieldnames = [('pscore_RAxML','pscoreout_RAxML'), \
				             ('pscore_RAxML_fixTreeEstBranch','fixTreeEstBranch.pscoreout_RAxML')]
		new_fieldnames = in_csv.fieldnames + [x[0] for x in pscore_fieldnames]
		out_csv = csv.DictWriter(out_f, new_fieldnames, delimiter=delimiter)
		# need to write out the field names
		#out_csv.writeheader()# lol this function only in 2.7 and i have 2.6 Orz 
		out_f.write(delimiter.join(new_fieldnames) + '\n')
		for obj in in_csv:
			motif_full_path = os.path.join(motif_dir, obj['motif'])
			print >> sys.stderr, motif_full_path
			for field, suffix in pscore_fieldnames:
				pscore_RAxML = os.popen("grep \"Total pair posterior\" {0}.{1}".format(motif_full_path, suffix)).read().strip()
				try:
					obj[field] = float( pscore_RAxML[len('Total pair posterior '):] )
				except ValueError:
					if len(os.popen("grep \"Fail to open\" {0}.{1}".format(motif_full_path, suffix)).read().strip()) > 0:
						obj[field] = 'NA'
			out_csv.writerow(obj)

if __name__ == "__main__":
	for filename in glob.iglob('motifs/*/*.fna.summary.pscore_added'):
		print >> sys.stderr, "annotating pscore_RAxML to {0}.RAxML_added....".format(filename)
		annotate_rank_summary_with_pscore_RAxML(filename)
