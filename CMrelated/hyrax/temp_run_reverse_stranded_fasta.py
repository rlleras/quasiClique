import os
import sys
import glob
from Bio import SeqIO
from miscParses import reverse_ID_strand

def main(dir_patterns):
	"""
	For every <subdir> in the directory,
	assume there is a <subdir>/<subdir>.fna
	Run reverse stranded on the fasta file
	Name the new fasta file XXX.fna.reversed
	"""
	for d in glob.iglob(dir_patterns):
		fna_filename = os.path.join(dir_name, d, d + '.fna')
		new_filename = fna_filename + '.reversed'
		with open(new_filename, 'w') as f:
		for r in SeqIO.parse(open(fna_filename), 'fasta'):
			f.write(">{id}\n{seq}\n".format(id=reverse_ID_strand(r.id), r.seq.reverse_complement()))
		os.system("mv {0} {0}.reversed".format(d))

if __name__ == "__main__":
	main(sys.argv[1])

