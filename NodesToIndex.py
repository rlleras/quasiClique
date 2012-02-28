from Bio import SeqIO # from BioPython package

class NodesToIndex:
	def __init__(self, d={}, highest_free_index=0, filename=None, igr=None):
		assert highest_free_index >= 0
		self.d = d
		self.highest_free_index = highest_free_index
		self.filename = filename
		self.igr      = igr
		if self.igr is not None:
			self.read_from_igr()

	def __getitem__(self, k):
		try:
			return self.d[k]
		except:
			self.d[k] = self.highest_free_index
			self.highest_free_index += 1
			return self.d[k]

	def clear(self):
		self.d = {}
		self.highest_free_index = 0
		self.igr = None
		self.filename = None

	def load(self, filename=None):
		"""
		load from <filename> or self.filename (if <filename> is None)
		"""
		assert filename is not None or self.filename is not None
		self.clear()
		if filename is None:
			self.filename = filename

		with open(self.filename) as f:
			self.igr = f.readline().strip() # first line in the igr filename
			for line in f:
				index, id = line.strip().split('\t')
				index     = int(index)
				assert index >= 0
				self.d[id] = index
				self.highest_free_index = max(self.highest_free_index, index)

	def write(self, filename=None):
		"""
		write out to <filename> or self.filename (if <filename> is None)
		"""
		assert filename is not None or self.filename is not None
		if filename is None:
			self.filename = filename

		with open(self.filename, 'w') as f:
			f.write( self.igr + '\n' )
			for id,index in self.d.iteritems():
				f.write("{index}\t{id}\n".format(index,id))

	def read_from_igr(self, igr):
		"""
		<igr> should be in FASTA format
		"""
		assert igr is not None or self.igr is not None
		self.clear()
		if igr is not None:
			self.igr = igr

		with open(self.igr) as f:
			for r in SeqIO.parse(f, 'fasta'):
				self.__getitem__(r.id)

