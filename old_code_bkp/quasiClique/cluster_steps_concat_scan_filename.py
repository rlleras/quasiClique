import os,re,sys
from cPickle import *
from Bio import SeqIO
from ushuffle import *

def concat_scan_filename(scan_filename, output_prefix, shuffled_copy):
	"""
	Assuming that the fasta file is used for CM scanning, and
	seq ids are in format <rb_fam>___<rb_id>. 

	Produces two files: <output_prefix>.as1chromosome.fna and
	<output_prefix>.as1chromosome.mapping_pickle, where the mapping
	pickle is a sorted list of (0-based start,0-based end,id)
	"""
	f = open(output_prefix+'.as1chromosome.fna', 'w')
	f.write(">FINAL\n")

	mapping = []

	cur_index = 0
	with open(scan_filename) as handle:
		for r in SeqIO.parse(handle, 'fasta'):
			seq_len = len(r.seq)
			seq = r.seq.tostring()
			f.write(seq)
			mapping.append( (cur_index, cur_index+seq_len-1, r.id) )
			shuffle1(seq, seq_len, 2)
			for blah in xrange(shuffled_copy):
				f.write(shuffle2())
			cur_index += (shuffled_copy + 1) * seq_len
	
	f.write('\n')
	with open(output_prefix+'.as1chromosome.mapping_pickle', 'w') as f:
		dump(mapping, f)

if __name__ == "__main__":
	concat_scan_filename(sys.argv[1], sys.argv[2], int(sys.argv[3]))
