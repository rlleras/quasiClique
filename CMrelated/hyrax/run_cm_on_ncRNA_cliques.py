import os,sys,glob
from cPickle import *

"""
 There should be a description.pickle that is a <motif_dir> --> dict info
 so we can get the family of each motif dir

 I've already created a 100X di-shuffled scan DB using Zasha's 2009-08-19-sto
 files for CM scanning stored universally on hyrax at /external3/home/etseng/CMscanDB_100X/
 the file names should be <family>.plus100Xdishuffled.fna

2010/08/16
 Goes through every motif file in every directory (look for the '.cm' files), 
 print the CM scan command according to the family this motif is. 
 The CM scan command runs only on top strand.
 The output is always the motif name appended with '.cmsearched'. 
"""

with open('description.pickle') as f: dir_info = load(f)

def main(d):
	fam = dir_info[d]['family']
	db_filename = "/external3/home/etseng/CMscanDB_100X/" + str(fam) + '.plus100Xdishuffled.fna'

	if not os.path.exists(db_filename):
		print >> sys.stderr, "NEED FILE {0}!!!!!!".format(db_filename)

	for file in glob.iglob("{d}/*.motif.*h[1-9]_[1-9]".format(d=d)):
		print("srun -p pubnorm -N 1 cmsearch --fil-T-hmm 10 -T 10 --tabfile {out} {cm} {db} > /dev/null &".format(\
			cm=file+'.cm',\
			db=db_filename,\
			out=file+'.cmsearched'))

if __name__ == "__main__":
	for d in os.listdir(.'):
		main(d)
		
