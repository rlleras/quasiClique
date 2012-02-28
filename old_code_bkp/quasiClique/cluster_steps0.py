import os,re,sys
from Bio import SeqIO, Seq

"""
Stuff for pre-BLASTing such as softmask-ing fasta files given SDUST results
"""

def mask_fasta_given_SDUST_intervals(fasta_filename, interval_filename, output_filename, rest_all_upper=False):
	"""
	Turn all positions indicated by SDUST's interval output file to LOWER CASE.
	Note that SDUST's positions are 0-based.
	SDUST's interval output looks like:
	>
	12 - 18
	345 - 351
	>
	>
	79 - 87

	If <rest_all_upper> is True, then regardless of what the original sequences cases were, turn everything
	that's not indicated by SDUST into upper case.
	"""
	f_out = open(output_filename, 'w')

	with open(interval_filename) as f_dust:
		with open(fasta_filename) as f_fasta:
			it = SeqIO.parse(f_fasta, 'fasta')
			
			r = it.next()
			f_dust.readline() # must be >
			to_mask = []

			for line in f_dust:
				if line.startswith('>'):
					m_seq = r.seq.tomutable()
					if rest_all_upper:
						m_seq = str(m_seq).upper()
					for s,e in to_mask:
						m_seq[s : e+1] = str(m_seq[s : e+1]).lower()
					# write out the sequence	
					f_out.write(">{id}\n".format(id=r.id))
					f_out.write("{s}\n".format(s=m_seq))

					r = it.next()
					to_mask = []
				else:
					to_mask.append( map(int, line.strip().split(' - ')) )

	f_out.close()


if __name__ == "__main__":
	mask_fasta_given_SDUST_intervals(*sys.argv[1:])
