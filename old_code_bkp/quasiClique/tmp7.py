import os,re,sys
from networkx import *
from networkx.readwrite import *
from cPickle import *
import cluster_steps as c
from miscMySQL import *

dirname = 'output/output_blast_etc/ALLFirm_RefSeq25_m30s0/'

conn = get_conn_Firm()
cursor = conn.cursor()
nodes_to_index = c.get_nodes_to_index_from_db(cursor)

G = Graph()
while G.number_of_edges() < 2100000:
	cursor.execute("select filename from files_todo order by rand() limit 1")
	try:
		filename = cursor.fetchone()[0]
	except:
		print >> sys.stderr, "couldn't fetch a file to process....no more?"
		break
	cursor.execute("delete from files_todo where filename=\"{0}\"".format(filename))
	print >> sys.stderr, filename
	c.step1_process_blast(dirname+filename, 35, nodes_to_index, G)
