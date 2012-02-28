import os,re,sys

"""
CMfinder often fails on fasta files with large number of seqs (ex: 100)
so go through the directories and print out all of them
"""
for d in os.listdir('.'):
	 if os.path.isdir(d):
		 if len(os.popen("ls {0}/{0}*motif.*".format(d)).read().strip()) == 0:
			 print d
