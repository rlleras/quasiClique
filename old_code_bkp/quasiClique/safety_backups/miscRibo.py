import os,re,sys,MySQLdb
from collections import defaultdict

db = MySQLdb.connect('xmen.cs.washington.edu','root','fidelcastro', 'ncRNA', unix_socket='/homes/gws/lachesis/mysqld-bio.sock', port=7777)
cursor = db.cursor()
DEFAULT_TABLE = 'Zasha_20081002_plus20071102_curated'

def get_all_rbtypes(table):
	"""
	    Simple returns a list of distinct familys in the entire table
	    ex: ['SAM','TPP',....]
	"""
	mysql_query = "select distinct family from " + table
	cursor.execute(mysql_query)
	result = []
	for r in cursor.fetchall():
		result.append(r[0])
	return result

def get_all_rbtype_and_count(table=DEFAULT_TABLE):
	"""
		Simple returns a dictionary of family --> rb_count 
		in the entire table. Probably used for evaluating
		CM scans, which is done on a DB containing ALL
		riboswitches...
	"""
	mysql_query = "select family,count(*) from {0} group by family" % table
	cursor.execute(mysql_query)
	result = {}
	for r in cursor.fetchall():
		result[r[0]] = int(r[1])
	return result

def get_ribo_count_by_acclist(acclist,table):
	"""
	    This is mainly used for estimating how much riboswitch there are
	    per accession list (however notice it's per ACCESSION not IGR)
	    Returns a dict of key(ribo_type) --> count(# of ribo of type in this acclist)
	"""
	# make the acclist to a sql-style string
	acclist = str(acclist)[1:-1]
	mysql_query = "select family,count(*) from %s where accession in (%s) group by family" % (table,acclist)
	#print mysql_query
	cursor.execute(mysql_query)
	result = {}
	for type,count in cursor.fetchall():
		result[type] = count
	return result

def get_ribo2(acc,start,end,table=DEFAULT_TABLE):
	"""
		Similar to ribo1 but different in the sense that:
		(1) DOESNT care about strand
		(2) returns ALL ribo hits if there are multiple ones
		(3) ONLY returns hits that are GOOD (>=50bp or >=50%), so there no marking for 'family-'

		This is usually used for IGR ribo hit check
	"""
	mysql_query = "select id,family,calc_overlap(start,end,%s,%s) as o,(end-start+1) from %s where acc=\"%s\" having o>=15 order by o desc" % (start,end,table,acc)
	cursor.execute(mysql_query)
	results = []
	for result in cursor.fetchall():
		overlap = int(result[2])
		rb_len = int(result[3])
		if overlap>=50 or overlap*1./rb_len>=0.5:
			results.append((result[0],result[1]))
	return results

def get_ribo1(acc,start,end,table=DEFAULT_TABLE):
	"""
		Query <table> to find matching ribos. Criteria for matching is
		" hits at least 50% or 50bp of the riboswitch "
		In case of hitting multiple riboswitches, return only 1 of them
		in the form of (rb_id,family).
		If it hits a ribo but with BAD overlap (failing the criteria), return (rb_id,family+'-')
	"""
	mysql_query = "select id,family,calc_overlap(start,end,%s,%s) as o,(end-start+1) as l from %s where acc=\"%s\" having o >= 15 order by o desc" % (start,end,table,acc);
	cursor.execute(mysql_query)
	result = cursor.fetchone()
	if result is None: return (None,None)
	else:
		overlap = int(result[2])
		rb_len  = int(result[3])
		if overlap>=50 or overlap*1./rb_len>=0.5:
			return (result[0],result[1])
		elif overlap>=15:
			return (result[0],result[1]+'-')
		else:
			return (None,None)

def get_phylum_of_ribo(rb_id,table):
	mysql_query = "select tax_phylum from Organism.organism where accession=(select accession from Ribo.%s where id=%s)" % (table,rb_id)
	cursor.execute(mysql_query)
	return cursor.fetchone()[0]

def read_phylum_distribution(filename='/homes/gws/lachesis/xmen/afterQual_pipeline/PythonScripts/BLASTparamTest/ncRNA_to_IGR_ratio/output/RefSeq25_IGRs_m30s0.ncRNA_phyla_distribution.txt'):
	"""
	Reads the distribution of ribo families per phylum file,
	returns a dictionary of {phylum} --> {family --> count}
	"""
	from collections import defaultdict
	d = defaultdict(lambda: {})
	f = open(filename, 'r')
	phyla = f.readline().strip().split('\t')
	for line in f:
		line = line.strip().split('\t')
		fam = line[0]
		for i,count in enumerate(line[1:]):
			if count == '': # means a count of 0
				d[phyla[i]][fam] = 0
			else:
				d[phyla[i]][fam] = int(count)
	f.close()
	return d

def get_fam_count(fasta_filename):
	from Bio import SeqIO
	from miscParses import parsed_accID
	total_fam_counts = defaultdict(lambda: 0)
        # read the scanned fasta so we know how many ribos per family there are
	# doing this everytime may be a waste, but it ensures we're having the right counts...
	for r in SeqIO.parse(open(fasta_filename), 'fasta'):
		print >> sys.stderr, "fasta reading....", r.id
		(acc, junk_version),strand,start,end = parsed_accID(r.id, True)
		(rb_id, rb_fam) = get_ribo1(acc, start, end)
                total_fam_counts[rb_fam] += 1
	return total_fam_counts
#print get_ribo1('JCVI_SCAF_1101667175582',225,275,1,"riboswitch_everything")
#print get_ribo1('NC_000964',1375612,1375799,-1,'riboswitch_everything')
