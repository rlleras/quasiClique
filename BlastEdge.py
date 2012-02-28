from interval import * # interval package

class BlastEdge:
	def __init__(self, i1, i2, intvl1, intvl2, score, opposite_strand):
		self.lines = {i1: IntervalSet([Interval(*intvl1)]), \
			      i2: IntervalSet([Interval(*intvl2)])}
		self.score = score
		self.muddy = False
		self.opposite_strand = opposite_strand
	
	def __str__(self):
		it = self.lines.iteritems()
		k1,v1 = it.next()
		k2,v2 = it.next()
		return "{0}:{1}\t{2}:{3}\t{4}({5})\t{6}".format(k1, v1, k2, v2, \
				self.score, self.muddy, self.opposite_strand)

	def lines_of(self, i):	
		return self.lines[i]
		
	def add(self, i1, i2, r1, r2, score):
		self.lines[i1].add(Interval(*r1))
		self.lines[i2].add(Interval(*r2))
		if score > self.score:
			self.score = score
			self.muddy = True

