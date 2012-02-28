from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Nexus.Trees import Tree
"""
	There should be supplements to the Biopython modules
"""
class FastaReader:
	"""
		This is meant to substitute for the Bio.SeqIO.to_dict method since some fasta files
		are too big to fit entirely to memory. The only requirement is that every id line
		begins with the symbol >. It is ok for the sequences to stretch multiple lines.
		The sequences, when read, are returned as Bio.SeqRecord objects.

		Example:
			r = FastaReader('output/test.fna')
			r['6C_49273_NC_008578/2259031-2259297'] ==> this shows the SeqRecord
	"""
	def __init__(self, fasta_filename):
		self.f = open(fasta_filename)
		self.d = {}

		while 1:
			line = self.f.readline()
			if len(line) == 0: break
			if line.startswith('>'):
				id = line.strip()[1:] # the header MUST be just 1 line
				#if id in self.d:
				#	print "duplicate id {0}!!".format(id)
				self.d[id] = self.f.tell()
				#print id,self.d[id]

	def __getitem__(self, k):
		if k not in self.d:
			raise Exception, "key {0} not in dictionary!".format(k)
		self.f.seek(self.d[k])
		content = ''
		for line in self.f:
			if line.startswith('>'):
				break
			content += line.strip()
		return SeqRecord(Seq(content), id=k)

	def keys(self):
		return self.d


class NewickReader:
	"""
		Just a wrapper around Bio.Nexus.Trees to read newick files. 
		In addition, since many of my newick taxon labels are just ncRNA <db_id>s, 
		support database lookup for these IDs as well.
	"""
	def __init__(self, filename):
		self.filename = filename
		self.tree = None

		f = open(self.filename,'r')
		chunk = f.read()
		f.close()
		self.tree = Tree(chunk)

	def distance(self, taxon1, taxon2):
		"""
			Note that here "taxon" simply means whatever the terminal nodes' data are.
			Since most of my newick files are labeled with <db_id>, it could just be ex: '34969'.
		"""
		id1 = self.tree.search_taxon(taxon1)
		id2 = self.tree.search_taxon(taxon2)
		return self.tree.distance(id1, id2)

def GCcontent(fasta_filename):
	acc,gc = 0,0
	for r in SeqIO.parse(open(fasta_filename), 'fasta'):
		s = r.seq
		acc += len(s)
		gc += s.count('G') + s.count('C')
	return gc*1./acc
