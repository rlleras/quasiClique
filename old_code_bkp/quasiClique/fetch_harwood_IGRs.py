import os,re,sys
import urllib2
import TableParse
from Bio import SeqIO, Seq

baseurl='http://www.pseudomonas.com/'

tax_to_acc = {\
		('aeruginosa', 'PAO1'): 'NC_002516',\
		('aeruginosa', 'PA7'):  'NC_009656',\
		('aeruginosa', 'LESB58'):'NC_011770',\
		('aeruginosa', 'PA14'): 'NC_008462',\
		('entomophila', 'L48'): 'NC_008027',\
		('fluorescens', 'PfO-1'):'NC_007492',\
		('fluorescens', 'Pf-5'): 'NC_004129',\
		('mendocina', 'ymp'):    'NC_009439',\
		('putida', 'W619'):      'NC_010501',\
		('putida', 'F1'):        'NC_009512',\
		('putida', 'KT2440'):    'NC_002947',\
		('stutzeri', 'A1501'):   'NC_009434',\
		('syringae', '1448A'):   'NC_005773',\
		('syringae', 'B728a'):   'NC_007005',\
		('syringae', 'DC3000'):  'NC_004578',\
}

def get_POG_link(locusID):
	url = baseurl + 'getAnnotation.do?locusID=' + locusID
	x = urllib2.urlopen( url ).read()
	i = x.find('searchAnnotation.do?')
	j = x.find('">', i)
	link = x[i:j].replace(' ','')

	i2 = x.find('POG', i, j)
	j2 = x.find('&operator', i, j)
	name = x[i2:j2]
	return name,link

def get_Upstream500Seq(locusID):
	url = baseurl + 'getAnnotation.do?locusID=' + locusID
	for x in TableParse.parse( urllib2.urlopen( url ).read() ):
		if len(x) > 1 and x[0]=='Upstream 500 BP Region':
			return x[1][ : x[1].find('\n')].replace(' ','')

genomic_location_rex = re.compile('(\d+)-(\d+)\((\S+)\)')
def make_ID_for_UPstream500Seq(locusID, tax_to_acc):
	url = baseurl + 'getAnnotation.do?locusID=' + locusID
	for x in TableParse.parse( urllib2.urlopen( url ).read() ):
		if len(x) > 1:
			if x[0]=='Strain':
				strain = x[1]
				if strain in tax_to_acc:
					acc = tax_to_acc[strain]
				else:
					acc = raw_input("Give the accession for {0}:".format(strain))
					tax_to_acc[strain] = acc
					print >> sys.stderr, "you typed:", acc
			elif x[0]=='Genomic location':
				m = genomic_location_rex.match( x[1].replace(' ','') )
				start,end,strand = int(m.group(1)), int(m.group(2)), m.group(3)
				print >> sys.stderr, "strand is ", strand
				if strand.find('+') >= 0:
					loc = "{0}/{1}-{2}".format(acc, start-500, start-1)
				else:
					loc = "{0}/{1}-{2}".format(acc, end+499, end+1)
				return loc

def getUpstream500_from_intergenic(intergenic_filename, protein_name, strand, start, end):
	"""
	Given <intergenic_filename> which is probably downloaded from Pseudomonas.com
	and is a fasta file of id format >{previous_protein}, {start}-{end}, {next_protein}
	
	If protein <strand> is '-', then we have to take the downstream 500bp
	(since it's actually upstream on - strand)

	Returns a Seq object containing the upstream seq,
	if return is None something	wacky happened!!!
	"""
	with open(intergenic_filename) as handle:
		records = list( SeqIO.parse(handle, 'fasta') )
		# because the id format is unusual, use the 'description' field for id-ing
		for i,r in enumerate(records):
			prev,loc,next = r.description.split(', ')
			_s, _e = map(int, loc.split('-'))
			if next == protein_name:
				# be careful: if it's - strand we don't want to include this seq
				if strand == '+':
					seq = r.seq.tostring()
					locs = [loc]
					while len(seq) < 500 and i > 0:
						i -= 1
						seq = records[i].seq.tostring() + seq
						locs.append( records[i].description.split(', ')[1] )
					seq = Seq.Seq( seq )
				else:
					seq = ''
					locs = []
					while len(seq) < 500 and i < len(records)-1:
						i += 1
						seq += records[i].seq.tostring()
						s,e = records[i].description.split(', ')[1].split('-')
						locs.append( "{0}-{1}".format(e,s) )
					# remember to reverse-complement it
					seq = Seq.Seq( seq ).reverse_complement()
					locs.reverse()
				return locs, seq
			elif prev == protein_name or _s > end:
				# this means we never saw a previous entry of xxx, xxx, <protein_name>
				# which may happen if <protein_name> overlaps with previous gene
				# be careful here: if it's + strand we don't want to include this seq
				if strand == '+':
					seq = ''
					locs = []
					while len(seq) < 500 and i > 0:
						i -= 1
						seq = records[i].seq.tostring() + seq
						locs.append( records[i].description.split(', ')[1] )
					seq = Seq.Seq( seq )
				else:
					seq = r.seq.tostring()
					locs = [loc]
					while len(seq) < 500 and i < len(records)-1:
						i += 1
						seq += records[i].seq.tostring()
						s,e = records[i].description.split(', ')[1].split('-')
						locs.append( "{0}-{1}".format(e,s) )
					seq = Seq.Seq( seq ).reverse_complement()
					locs.reverse()
				return locs, seq
	return None,None	

def get_POGs(link):
	POGs = []
	url = baseurl + link
	x = urllib2.urlopen( url ).read()
	for p in TableParse.parse( x ):
		if len(p) > 1 and p[0] == '-->':
			POGs.append( p[1] )
	return POGs

#def main(locusID, f):
#	POG_link = get_POG_link(locusID)
#	POGs = get_POGs( POG_link )
#	for pog in POGs:
#		seq = get_Upstream500Seq( pog )
#		f.write(">500BP_UPSTREAM_{0}\n{1}\n".format(pog, seq))

def main(locusID, f):
	POG_name, POG_link = get_POG_link(locusID)
	POGs = get_POGs( POG_link )
	for pog in POGs:
		url = baseurl + 'getAnnotation.do?locusID=' + pog
		for x in TableParse.parse( urllib2.urlopen( url ).read() ):
			if len(x) == 0: 
				continue
			if x[0] == 'Strain':
				strain = x[1][ : x[1].find('(')].strip().split()
				key = (strain[1],strain[2])
				if key in tax_to_acc:
					acc = tax_to_acc[key]
				else:
					acc = raw_input("Give the accession for {0}:".format(key))
					tax_to_acc[key] = acc
					print >> sys.stderr, "you typed:", acc
			elif x[0] == 'Genomic location':
				m = genomic_location_rex.match( x[1].replace(' ','') )
				start,end,strand = int(m.group(1)), int(m.group(2)), m.group(3)
				strand = '+' if strand.find('+') >= 0 else '-'
				break

		strain = 'P.' + strain[1] + strain[2]
		locs,seq = getUpstream500_from_intergenic(strain+'.intergenic', pog, strand, start, end)
		if seq is None:
			print >> sys.stderr, "something wacky happened to seq retrieving for {0}({1}) from file {2}!!!".format(\
					pog, strand, strain+'.intergenic')
		else:
			f.write(">" + strain + '_' + acc)
			f.write("_" + ",".join(locs))
			f.write("_" + POG_name + "_" + pog + "\n")
			f.write("{0}\n".format(seq))
			f.flush()

if __name__ == "__main__":
#	tax_to_acc_cache = {}
#	for d in os.listdir('motif_storage/'):
#		print >> sys.stderr, d
#		dd = os.path.join( 'motif_storage/', d)
#		for r in SeqIO.parse(open(dd+'/'+d+'.fna'), 'fasta'):
#			if not r.id.startswith('500BP_UPSTREAM_'): continue
#			id = r.id[len('500BP_UPSTREAM_'):]
#			id2 = make_ID_for_UPstream500Seq(id, tax_to_acc_cache)
#			id2 = id2.replace("/", "\\/")
#			print >> sys.stderr, "new id for {0} is {1}".format(id, id2)
#			os.system("sed -i 's/500BP_UPSTREAM_{0}/{1}/g' {2}/*.*".format(id, id2, dd))

	output_filename = sys.argv[1]
	handle = open(output_filename, 'w')
	with open('../bioinfo/harwood_Pseudomonas_aeruginosa_cdiGMPexpressed.justLocusTag.txt') as f:
		for line in f:
			id = line.strip()
			main(id, handle)
	handle.close()
