import os,re,sys,fnmatch,random
from Bio import SeqIO
from ushuffle import *
from time import time

random.seed(time())

d = '/homes/gws/lachesis/Larry_Riboswitch/Latest_ribo_20081002_plus20071102_curated/fnas'
select_ratio = .2

# first we write a FASTA file that contains <select_ratio> percent of ALL ncRNAs, 
# shuffled while perserving dinucleotide frequency
line_count = 0
with open('./output/random_selected_ncRNA_dishuffled.fna','w') as outhand:
	for filename in fnmatch.filter(os.listdir(d),'*.fna'):
		print >> sys.stderr, filename+'........randomly selecting then dishuffling....'
		with open(d+'/'+filename) as handle:
			for r in SeqIO.parse(handle,'fasta'):
				if random.random() > select_ratio:
					continue
				outhand.write('>'+r.id+'\n')
				s = r.seq.tostring()
				shuffle1(s,len(s),2)	
				for x in xrange(random.randint(1,10)):
					s1 = shuffle2(s,len(s),2)
				outhand.write(s1+'\n')
				line_count += 1

# run WU-BLAST on this dishuffled FASTA
os.chdir('./output/')
os.system("xdformat -n -o {0} {1}".format('random_selected_ncRNA_dishuffled.fna','random_selected_ncRNA_dishuffled.fna'))
for i in xrange(line_count-1):
	os.system("blastn random_selected_ncRNA_dishuffled.fna random_selected_ncRNA_dishuffled.fna -M 5 -N -4 -Q 10 -R 6 -W 7 -E 10 -noseqs -cpus 4 -mformat 2 -o random_selected_ncRNA_dishuffled.M5N4Q10R6W7E10.WUblast_part{0} -qrecmin {1} -qrecmax {2} -dbrecmin {3}".format(i,i+1,i+1,i+2))
os.chdir('../')
