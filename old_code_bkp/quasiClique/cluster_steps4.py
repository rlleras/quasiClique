import os,re,sys,fnmatch,operator
from cPickle import *
from cluster_steps_settings import *
from Bio import SeqIO
from miscParses import parsed_accID
from miscRibo import get_ribo1
from miscUtils import ToDirAndBack
from collections import defaultdict
from operator import itemgetter
import miscCMF

"""
Used to process motifs that have already been rm_dup and comb_motif-ed.
Depending on the case, the motif's should run through rm_dup_instance.pl
 first, with the resulting motifs having suffix .dup_rmed
Uses Zizhen's rank_cmfinder.pl script to rank motifs, then uses
 Zasha's scripts to unblock, color the motifs, then output to HTML format.
"""

def read_cmfinder_motif(motif_filename):
	"""
	Reads a motif file (outputted by CMfinder, stockholm format)
	Returns (ncRNA family, # of motif members belonging to the family, # of motif members)
	ncRNA family is determined by plurality. family is None if plurality are not ncRNAs.
	"""
	tally_by_family = defaultdict(lambda: 0)
	motif_size = 0
	with open(motif_filename) as f:
		f.readline()
		f.readline()
		for line in f:
			if len(line.strip()) == 0:
				continue
			feature, text, annot = line.strip().split(None,2)
			if feature == '#=GS' and annot.startswith('DE'):
				break
			if feature == '#=GS' and annot.startswith('WT'):
				motif_size += 1
				print >> sys.stderr, "looking up id", text
				#id = id[:id.rfind('_')] # what was this for?
				(acc,junk),strand,start,end = parsed_accID(text,True)
				ncRNA_id, ncRNA_family = get_ribo1(acc,start,end)
				tally_by_family[ncRNA_family] += 1
	print >> sys.stderr, "motif filename is", motif_filename
	print >> sys.stderr, "tally by family is", tally_by_family
	tally_by_family = tally_by_family.items()
	tally_by_family.sort(key=itemgetter(1))
	fam,count = tally_by_family[-1]
	return fam,count,motif_size

def eval_original_fna(fasta_filename):
	"""
	Reads a fasta file and returns (ncRNA family, # of seqs belonging to the family, clique size)
	ncRNA family is determined by plurality, which can be None.
	"""
	tally_by_family = defaultdict(lambda: 0)
	ids_hit_by_family = defaultdict(lambda: set())
	clique_size = 0
	for id in os.popen("grep \"^>\" " + fasta_filename):
		id = id.strip()[1:]
		clique_size += 1
		#id = id[:id.rfind('_')] # what was this for???
		(acc,junk),strand,start,end = parsed_accID(id,True)
		ncRNA_id, ncRNA_family = get_ribo1(acc,start,end)
		tally_by_family[ncRNA_family] += 1
		ids_hit_by_family[ncRNA_family].add( ncRNA_id )
	tally_by_family = tally_by_family.items()
	tally_by_family.sort(key=itemgetter(1))
	fam,count = tally_by_family[-1]
	# HACK HERE!!!
	if fam is None and len(tally_by_family) > 1:
		if tally_by_family[-2][1] >= 0.5*clique_size:
			lesser_fam = tally_by_family[-2][0] + '-'
			return lesser_fam,tally_by_family[-2][1],clique_size,ids_hit_by_family[tally_by_family[-2][0]]
		elif tally_by_family[-2][1] >= 3:
			lesser_fam = tally_by_family[-2][0] + '--'
			return lesser_fam,tally_by_family[-2][1],clique_size,ids_hit_by_family[tally_by_family[-2][0]]
	return fam,count,clique_size,ids_hit_by_family[fam]

def eval_clique(Q, cursor):
	"""
	Given Q which is a clique containing node indices, look it up
	on the db using cursor.
	
	Like eval_original_fna, returns:
	(fam, # of seqs belonging to fam, clique_size, fam ids hit)
	"""
	FETCH_SQL = "SELECT n.id,s.start,s.end \
				FROM sets_for_nodes s \
				LEFT JOIN nodes_to_index AS n \
				ON (s.nodes_ind=n.ind) WHERE i={i}"

	T = defaultdict(lambda: 0)       # fam ---> hit count
	H = defaultdict(lambda: set())   # fam ---> set of ids hit

	for i in Q:
		cursor.execute( FETCH_SQL.format(i=i) )
		_id,_loc_start,_loc_end = cursor.fetchone()
		(acc,junk),strand,start,end = parsed_accID(_id, True, _loc_start, _loc_end)
		id, fam = get_ribo1( acc, start, end )
		T[fam] += 1
		H[fam].add( id )

	T = T.items()
	T.sort(key=itemgetter(1), reverse=True)
	fam,count = T[0]
	# HACK
	if fam is None and len(T)>1:
		if T[1][1] >= 0.5*len(Q):
			return T[1][0]+'-',T[1][1],len(Q),H[T[1][0]]
		elif T[1][1] >= 3:
			return T[1][0]+'--',T[1][1],len(Q),H[T[1][0]]
	return fam, count, len(Q), H[fam]

class MotifRankInfo:
	NO_MOTIF_RANK = -1

	def __init__(self, fam, count, motif_size, motif_filename, fam_o, count_o, clique_size, rank):
		self.fam = fam
		self.count = count
		self.motif_size = motif_size
		self.motif_filename = motif_filename
		self.fam_o = fam_o
		self.count_o = count_o
		self.clique_size = clique_size
		self.rank = rank
		self.overwrite_bgcolor = None

	def __str__(self):
		return """
		motif rank  : {rank}
		motif file  : {motif_filename}
		motif family: {fam}
		motif spec  : {count}/{motif_size}
		clique family: {fam_o}
		clique spec  : {count_o}/{clique_size}
		""".format(\
				rank=self.rank,\
				motif_filename=self.motif_filename,\
				fam=self.fam,\
				count=self.count,\
				motif_size=self.motif_size,\
				fam_o=self.fam_o,\
				count_o=self.count_o,\
				clique_size=self.clique_size)

	def str_for_table(self, linkdir, bgcolor):
		if self.overwrite_bgcolor is not None:
			bgcolor = self.overwrite_bgcolor
		return """
		<tr bgcolor="{bgcolor}">
		  <td><a href="{link}">{ID}(<a href="{sto}">.sto</a>)</a></td>
		  <td>{rank}</td>
		  <td>{fam}</td>
		  <td>{count}/{motif_size}</td>
		  <td>Was a {count_o}/{clique_size} {fam_o} clique</td>
		</tr>
		""".format(bgcolor=bgcolor,\
				ID=os.path.basename(self.motif_filename),\
				link=os.path.join(linkdir,self.motif_filename),\
				sto=os.path.join(linkdir,self.motif_filename[:-len('.html')]+'.unblocked'),\
				rank=self.rank,\
				fam=self.fam,\
				count=self.count,\
				motif_size=self.motif_size,\
				count_o=self.count_o,
				clique_size=self.clique_size,\
				fam_o=self.fam_o)

	def table_header(self):
		return """
		<table border=1 cellspacing=4 cellpadding=4>
		<tr>
		  <td>ID</td>
		  <td>Rank</td>
		  <td>Family</td>
		  <td># of known ncRNAs in motif</td>
		  <td>Original clique was</td>
		</tr>
		"""

	def table_tail(self):
		return """
		</table>
		"""

def rank_motifs(dir_of_motif_dir, webdir, use_pscore=True, use_rankpl=False):
	"""
	Exactly one of use_pscore/use_rankpl should be True.
	"""
	assert os.path.isdir(webdir)
	if not operator.xor(use_pscore, use_rankpl):
		raise Exception, "exactly use_pscore/use_rankpl should be True!"

	total_ranks = []
	
	dir_of_motif_dir = os.path.abspath( dir_of_motif_dir )

	for d in os.listdir(dir_of_motif_dir):
		dd = os.path.join( dir_of_motif_dir, d )
		if not os.path.isdir(dd):
			continue
		print >> sys.stderr, "ranking motifs for directory {0}...".format(d)

		with ToDirAndBack(dd):
			fam_o,count_o,clique_size,ids_hit = eval_original_fna( d + '.fna' )

			if use_rankpl:
				ranks = miscCMF.rank_cmfinder_score( dd )
			else: # use pscore
				ranks = miscCMF.pfold_pscore( dd )
		
			if len(ranks) == 0:
				if fam_o is not None:
					# if this clique originally was a ncRNA clique, if so, 
					# we put it in total_ranks but with a rank of -1.0
					total_ranks.append( MotifRankInfo(fam='NA',\
							count=0,\
							motif_size=0,\
							motif_filename=d+'.fna',\
							fam_o=fam_o,\
							count_o=count_o,\
							clique_size=clique_size,\
							rank=MotifRankInfo.NO_MOTIF_RANK) )
					print >> sys.stderr, "copying fastafile {0} to webdir {1}".format(d+'.fna',webdir)
					os.system("cp {fasta} {webdir}".format(fasta=d+'.fna', webdir=webdir))
				continue

			rank_index,motif_filename = ranks[0] # for now take just the top rank motif
			fam,count,motif_size = read_cmfinder_motif(motif_filename)			
			total_ranks.append( MotifRankInfo(fam=fam,\
					count=count,\
					motif_size=motif_size,\
					motif_filename=os.path.basename(motif_filename)+'.html',\
					fam_o=fam_o,\
					count_o=count_o,\
					clique_size=clique_size,\
					rank=rank_index) )

			# use Zasha's script to create a colorful alignment html, then copy it to web directory
			os.system("perl ~/lib/perl/StockholmUnblock.pl {0} tmp.unblocked".format(motif_filename))
			os.system("cmzasha --GSC-weighted-consensus tmp.unblocked tmp 3 0.97 0.9 0.75 4 0.97 0.9 0.75 0.5 0.05")
			os.system("perl -I/homes/gws/lachesis/lib/perl ~/lib/perl/FancifyStockholm.pl "\
					"-noWarnParseHitId -forNewMotifsWeb -noURL -highlightCovarying tmp {0}.html".format(motif_filename))
			os.system("cp {motif}.html {webdir}".format(motif=motif_filename, webdir=webdir))
			os.system("cp tmp.unblocked {webdir}/{motif}.unblocked".format(motif=os.path.basename(motif_filename),webdir=webdir))
		
	return total_ranks

def color_rank_in_html(total_ranks, html_filename, linkdir):
	"""
	Outputs to <output> an HTML table that shows the rankings for motifs, along with 
	other info generated by Zizhen's rank_cmfinder.pl script.
	"""
	total_ranks.sort(key=lambda obj:obj.rank, reverse=True)

	with open(html_filename, 'w') as f:
		f.write("""
		<html>
		<body>
		""")

		# write the header
		f.write( total_ranks[0].table_header() + '\n' )

		for rank_obj in total_ranks:
			if rank_obj.fam is None:
				if rank_obj.fam_o is None:
					bgcolor = '#FFFFFF'
				else:
					bgcolor = '#FF66CC'
			else:
				bgcolor = '#99CCFF'
			f.write( rank_obj.str_for_table(linkdir, bgcolor) + '\n' )

		f.write( total_ranks[0].table_tail() + '\n' )

		f.write("""
		</body>
		</html>
		""")

def compare_ranking_pickles_usingSpearman(pickle_filename1, pickle_filename2):
	"""
	Compare two ranking pickles (produced by rank_motifs) which is a 
	(not always sorted) list of MotifRankInfo using Spearman correlation. 

	Assumes MotifRankInfo obj's motif_filename is of form:
	    grid05_168_173.fna_split03.1.motif.h1_3.h2_4.h2_3.h2_5.dup_rmed.html

	So the key we will be grabbing is grid05_168_173
	"""
	from itertools import groupby
	from math import sqrt
	# inner function for taking a list of MotifRankInfo objs and
	# returning a dictionary of mapping obj.motif_filename --> rank
	# in cases of tie, must assign a rank that's the avg of their ranks
	# ex: 1,2,2,4 --> 1,2.5,2.5,4
	def assign_rank(MRIs):
		result = {}
		MRIs.sort(key=lambda obj:round(obj.rank), reverse=True)
		cur_rank = 1
		for _duncare,_iter in groupby(MRIs, lambda obj:obj.rank):
			to_assign = [obj.motif_filename[:obj.motif_filename.find('.fna')] for obj in _iter]
#			print >> sys.stderr, "for obj.rank={0} there are {1}".format(_duncare, to_assign)
			# so the rank for these equal-valued to_assign is [cur_rank .. cur_rank+len(to_assign)-1]
			_rank_to_assign = (cur_rank + cur_rank + len(to_assign) - 1) / 2.
#			print >> sys.stderr, "assigning them rank {0}....".format(_rank_to_assign)
			for k in to_assign:
				# hack for part2/part1
				k = k.replace('part2','part1')
				k = k[:k.find('_size')]
				result[ k ] = _rank_to_assign
			cur_rank += len(to_assign)
#			raw_input()
		return result

	with open(pickle_filename1) as handle:
		A1 = assign_rank( load(handle) )
	with open(pickle_filename2) as handle:
		A2 = assign_rank( load(handle) )
	
	sum_of_xy = 0
	sum_of_x  = 0
	sum_of_x_squared = 0
	sum_of_y = 0
	sum_of_y_squared = 0
	n = 0
	# use Pearson's correlation coefficient with ranks
	omitted = 0
	for k in A1:
		if k not in A2:
			omitted += 1
			print >> sys.stderr, "omitting key {0} because not in A2".format(k)
		else:
			n += 1
			x, y = A1[k], A2[k]
			sum_of_xy += x*y
			sum_of_x  += x
			sum_of_x_squared += x**2
			sum_of_y  += y
			sum_of_y_squared += y**2

	rho = n*sum_of_xy - sum_of_x*sum_of_y
	rho /= sqrt(n*sum_of_x_squared - sum_of_x**2) * sqrt(n*sum_of_y_squared - sum_of_y**2)

	print >> sys.stderr, "{0} omitted".format(omitted)
	return rho

if __name__=="__main__":
#	print compare_ranking_pickles_usingSpearman('split1_rankpl.pickle','split2_rankpl.pickle')
#	sys.exit(-1)
	from optparse import OptionParser
	import platform
	parser = OptionParser()
	parser.add_option("-d", "--dir", help="Directory that contains the directories of motifs")
	parser.add_option("-w", "--webdir", help="Web directory that motifs will be copied to, ex: ~/www/quasiClique/Firm/motifs/")
	parser.add_option("-l", "--linkdir", help="Link directory that's used for formatting hyperlinks, ex: http://www/homes/lachesis/quasiClique/Firm/motifs/")
	parser.add_option("-r", "--rankmethod", help="Method for ranking, options are pscore[def] and rankpl", default="pscore")
	parser.add_option("-p", "--pickle", help="Output pickle filename for total_ranks")
	parser.add_option("-o", "--html", help="Output html filename")
	(options, args) = parser.parse_args()

	print("CONNECTION: {conn}".format(conn=CONN_FUNC))
	if os.path.exists(options.pickle):
		print("Output pickle file already exists! This will overwrite it!")
	if os.path.exists(options.html):
		print("Output html file already exists! This will overwrite it!")
	raw_input("--- is this OK?")

	if options.rankmethod == 'pscore':
		use_pscore,use_rankpl = True,False
	elif options.rankmethod == 'rankpl':
		use_pscore,use_rankpl = False,True
	else:
		raise Exception, "Invalid rankmethod choice!"

	total_ranks = rank_motifs(options.dir, options.webdir, use_pscore, use_rankpl)

	with open(options.pickle, 'w') as f:
		dump(total_ranks,f)

	color_rank_in_html(total_ranks, options.html, options.linkdir)

