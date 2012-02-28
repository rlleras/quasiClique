import os,re,sys,fnmatch
from ushuffle import shuffle
from Bio import SeqIO
from urllib import urlretrieve
from miscParses import parsed_accID
from miscMySQL import *
from cPickle import *


EUTIL_URL = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=genome&retmode=text&rettype=fasta&'

def padseq(input_filename,output_dirname,tablename='Zasha_20081002_plus20071102_curated'):
	"""
		If we can get the real flank by querying NCBI, then get it.
		Else, pad with shuffled version (preserving dinucleotide freq).
	"""
	conn = get_conn_ncRNA()
	cursor = conn.cursor()

	output_prefix = output_dirname+'/'+os.path.basename(input_filename)[:-4]
	f = open(output_prefix+'_padseqlen.fna','w')
	map_padded_to_orig_id = {}

#	rex = re.compile('(\S+)\\/(\d+)-(\d+)')
	tmp_outfile = os.tempnam()
	for r in SeqIO.parse(open(input_filename),'fasta'):
		i = r.id.rfind('_')
		family,db_id = r.id[:i],r.id[(i+1):]
		cursor.execute("select acc,start,end,strand from {0} where id={1}".format(tablename,db_id))
		try:
			acc,start,end,strand = cursor.fetchone()
		except:
			print("failed on {0}".format(db_id))
			continue
#		print >> sys.stderr, db_id,acc,start,end,strand
#		(acc,duncare),start,end,strand = parsed_accID(r.id)
		seqlen = end-start+1
		final_start = max(1,start-seqlen)
		final_end = end+seqlen

		map_padded_to_orig_id["{0}/{1}-{2}".format(acc,final_start,final_end)] = r.id

		url = EUTIL_URL + 'id=' + acc + '&strand=' + str(strand) + '&seq_start=' + str(final_start) + '&seq_stop=' + str(final_end)
		urlretrieve(url, tmp_outfile)
	
		try:
			rec = SeqIO.parse(open(tmp_outfile),'fasta').next().seq.tostring()
		except:
			rec = None
		if rec is None or len(rec)<>(final_end-final_start+1):
			s = r.seq.tostring()
			tmp_front = ''
			tmp_back = ''
			tmp_front = shuffle(s,seqlen,2)
			tmp_back  = shuffle(s,seqlen,2)
			f.write(">{4}_{0}/{1}-{2}\n{3}\n".format(acc,final_start,final_end,tmp_front+s+tmp_back,r.id))
		else:
			f.write(">{4}_{0}/{1}-{2}\n{3}\n".format(acc,final_start,final_end,rec,r.id))
	f.close()

	f = open(output_prefix+'_padseqlen.map','w')
	dump(map_padded_to_orig_id,f)
	f.close()

def padseq_by_directory(input_dirname,output_dirname):
	for filename in fnmatch.filter(os.listdir(input_dirname),'*.fna'):
		padseq(input_dirname+'/'+filename,output_dirname)

if __name__=="__main__":
	padseq_by_directory(sys.argv[1],sys.argv[2])
