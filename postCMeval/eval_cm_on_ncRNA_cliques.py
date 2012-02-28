from settings import *
from miscInfernal import read_tab_cmscanned

"""
 Needs: description.pickle, 
        /external3/home/etseng/CMscanDB_100X/, (or any pre-set dir)
        run_cm_on_ncRNA_cliques.py

 Reads the .cmsearched files and compiles them into a single print-out
 that is tab-delimited by:
	family, family size, motif name, rank score, pscore, TP count,...
        TP scores(comma-delimited), FP scores(comma-delimited)

 FPs are always identified with the FAKE_ seq id prefix!

2010/08/16
 This was run after {run_cm_on_ncRNA_cliques.py} to see the correlation
 between pscore/rank scores and CM scan results. 
"""

def main():
	os.chdir(MOTIF_DIR)
	fieldnames = ['FAMILY', 'FAMSIZE', 'MOTIF', 'RANK', 'PSCORE', 'TOTAL_RECALL', 'TPS', 'FPS']
	out_csv = csv.writer(sys.stdout, fieldnames, delimiter='\t') # NOTE THAT we CAN't USE comma as delim becuz we use then for TP,FP scores!
	out_csv.writerow(fieldnames)
	for d in os.listdir(os.path.curdir):
		fam = dir_info[d]['family']
		motif_info = {}
		# read the .fna.summary.pscore_added file
		for obj in csv.DictReader(\
				open(os.path.join(d, d + '.fna.summary.pscore_added')),\
				delimiter=','):
			motif_name = os.path.basename(obj['motif']) 
			motif_info[motif_name] = obj
		for file in glob.iglob(d + "/*.cmsearched"):
			print >> sys.stderr, "reading {0}....".format(file)
			tp_count = 0
			tp = []
			fp = []
			for id,(score,start,end) in read_tab_cmscanned(file, 'score', 0).iteritems():
				if id.startswith('FAKE_'): # is a FP!
					fp.append( str(score) )
				else:
					tp.append( str(score) )
					tp_count += 1
			motif_name = os.path.splitext(os.path.basename(file))[0]
			out_csv.writerow([fam, \
					db_summary[fam]['TRUE'], \
					motif_name, \
					motif_info[motif_name]['Rank.index'],\
					motif_info[motif_name]['pscore'],\
					tp_count,\
					",".join(tp),\
					",".join(fp)])

if __name__ == "__main__":
	main()
