import os,re,sys

def squid_shuffle_seq(input_filename,output_filename,n=10,prefix='shuffled_'):
	"""
		Sean Eddy's squid shuffle program doesn't rename the seq IDs, 
		so here we rename them
	"""
	dummy_count = 1
	f = open(output_filename,'w')
	for line in os.popen("shuffle -d -n {0} {1}".format(n,input_filename)):
		if line.startswith('>'):
			f.write(">"+prefix+str(dummy_count)+"_"+line.strip()[1:]+'\n')
			dummy_count += 1
		else:
			f.write(line.strip()+'\n')
	f.close()

def uniform_distribute_fasta(fasta_filename,rex,output_filename,size_limit=float("Inf")):
	"""
		This is used to uniformly distribute the sequences from a fasta file...
	"""
	from Bio import SeqIO
	from collections import defaultdict
	binned = defaultdict(lambda: [])
	with open(fasta_filename) as f:
		r_dict = SeqIO.to_dict(SeqIO.parse(f,'fasta'))
		for k,v in r_dict.iteritems():
			m = rex.match(k)
			binned[m.group(1)].append(v)
	lst = binned.keys()
	not_empty = True
	f = open(output_filename,'w')
	while not_empty:
		not_empty = False
		for k in lst:
			if len(binned[k])!=0:
				not_empty = True
				x = binned[k].pop()
				s = x.seq.tostring()
				if len(s) > size_limit or len(s) < 150: continue
				f.write(">{0}\n{1}\n".format(x.id,s))
	f.close()
				
if __name__=="__main__":
#	squid_shuffle_seq(sys.argv[1],sys.argv[2])
#	uniform_distribute_fasta('ALLFirmicutes_20081002_plus20071102_curated_withDup_modishuffled.fna',re.compile('\S+_\d+_(\S+)_\d+'),'ALLFirmicutes_20081002_plus20071102_curated_withDup_modishuffled.uniform.fna')
	uniform_distribute_fasta('picked_random_nonribo.fna',re.compile('random_(\S+)\/\S+'),'picked_random_nonribo.seqlen150_1000.uniform.fna',1000)
