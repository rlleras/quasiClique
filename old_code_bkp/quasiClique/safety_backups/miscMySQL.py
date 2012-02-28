#import pymysql
#
#def get_conn_Actino():
#	conn = pymysql.connect(host='xmen', port=7777, user='root', password='fidelcastro', db='ALLActino')
#	return conn

import MySQLdb
from MySQLdb.cursors import DictCursor

ncRNA_tablename = 'Zasha_20081002_plus20071102_curated'

def get_conn_ncRNA():
	conn = MySQLdb.connect(host='xmen',port=7777,user='root',passwd='fidelcastro',db='ncRNA')
	return conn

def get_conn_NCBI():
	conn = MySQLdb.connect(host='xmen',port=7777,user='root',passwd='fidelcastro',db='NCBI')
	return conn

def get_conn_Actino():
	conn = MySQLdb.connect(host='xmen',port=7777,user='root',passwd='fidelcastro',db='ALLActino')
	return conn

def get_conn_Firm():
	return MySQLdb.connect(host='xmen',port=7777,user='root',passwd='fidelcastro',db='ALLFirm')

def get_conn_Upwelling():
	return MySQLdb.connect(host='xmen',port=7777,user='root',passwd='fidelcastro',db='Upwelling')

def get_dict_cursor(conn):
	return DictCursor(conn)

def get_families(phylum=None,tablename='Zasha_20081002_plus20071102_curated'):
	conn = get_conn_ncRNA()
	cursor = conn.cursor()

	if phylum is None:
		cursor.execute("select distinct family from {0}".tablename)
	else:
		cursor.execute("select distinct z.family as family from {0} as z left join NCBI.cache_acc_to_tax as c on (z.acc=c.acc) where c.phylum='{1}'".format(tablename,phylum))
	
	return map(lambda x: x[0],cursor.fetchall())


def gen_fasta_from_db(phylum):
	with get_conn_ncRNA() as conn:
		cursor = get_dict_cursor(conn)
		cursor.execute("SELECT z.family as family,z.id as id,strand,start,end,seq \
						FROM Zasha_20081002_plus20071102_curated as z \
						LEFT JOIN NCBI.cache_acc_to_tax as c \
						ON (z.acc=c.acc) WHERE phylum='{0}'".format(phylum))
		for r in cursor.fetchall():
			print(">{0}_{1}\n{2}".format(r['family'],r['id'],r['seq']))
	
