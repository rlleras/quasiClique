import os,sys,glob
import myUShuffle, miscInfernal
from Bio import SeqIO

convert_py = '/home/etseng/lib/python/convert/sto_to_fasta.py'
db_func = myUShuffle.create_scan_db
tally_func = miscInfernal.read_tab_cmscanned #(tab_filename, best_e_or_score, cutoff)
X = 100
"""
 Given a .sto file which is Zasha's ncRNA family stockholm file
 1) convert to .fna
 2) run CMfinder on it
 3) create a scan db with (X-1) times dishuffled negative controls
 4) run every .cm file on the scan db
 5) find the CM score (or E-value?) cutoff that gives max recall
    with FDR < 1e-3
"""

def main(sto_filename):
	fna_filename = os.path.splitext(sto_filename)[0] +'.fna'
	os.system("python " + convert_py + " " + sto_filename)
	os.system("cmfinder.pl -def " + fna_filename)
	db_filename = db_func(fna_filename, X)

	# run pscore-pfold
	os.system("bash run_pscore.sh")
	# run ranking score
	os.system("rank_cmfinder.pl -w -rank \"{fna}.[1-9].motif.*\" {fna}.summary".format(fna=fna_filename))
	# read the summary file
	rank_of = {} # motif file --> rank score (in string form)
	with open(fna_filename+'.summary') as h:
		h.readline() #header
		for line in h:
			x = line.strip().split(',')
			rank_of[ os.path.basename(x[0]) ] = x[-1]
	
	pscore_of = {} # motif file --> pscore (in string form)
	motif_files = glob.glob('*.fna.*.motif.*h[1-9]_[1-9]')
	# run cmbuild
	for file in motif_files: os.system("cmbuild {file}.cm {file}".format(file=file))
	# run cmscan
	f = open(fna_filename + '.CMscan.out', 'w')
	f.write("motif\trank\tpscore\tTPhits\tFPhits\n")
	for file in motif_files:
		pscore_of[file] = os.popen("grep \"Total pair posterior\" {file}.pscoreout".format(file=file)).read().strip()[len('Total pair posterior'):]
		out_filename = file+'.cmsearched'
		os.system("cmsearch --toponly --tabfile {out} {cm} {db}".format(\
				cm = file+'.cm',\
				db = db_filename,\
				out= out_filename))#
		tp_hits = []
		fp_hits = [] 
		for id,(score,start,end) in tally_func( out_filename, 'score', cutoff=0 ).iteritems():
			if id.startswith('FAKE_'): # is a false positive hit!
				fp_hits.append( str(score) )
			else:
				tp_hits.append( str(score) )
		# for CM file <file>, fp_hits and tp_hits each contain list of TP/FP hit scores
		f.write( file + '\t' + str(rank_of[file]) + '\t' + str(pscore_of[file]) + '\t' + ",".join(tp_hits) + '\t' + ",".join(fp_hits) + '\n')
	f.close()

	# now expand the .summary file so that it contains the extra pscore column
	f = open(os.tempnam(), 'w')
	h = open(fna_filename + '.summary')
	f.write( h.readline().strip() + ',pscore\n' )
	for line in h:
		motif_name = os.path.basename( line.split(',')[0] )
		f.write( line.strip() + ',' + pscore_of[motif_name] + '\n')
	h.close()
	f.close()
	os.system("mv {old} {new}".format(old=f.name, new=h.name))

	
if __name__ == "__main__":
	main( sys.argv[1] )
				

