import os,re,sys
from cluster_steps_settings import *
#from scipy import histogram
#import dislin
"""
	These are mainly for calculating statistics & evaluating results
	from the approximate-clique finding through MySQL.
	For now it assumes it works on the database 'ALLFirm'
"""
def histogram_of_node_lengths():
	conn = get_conn_Firm()
	cursor = conn.cursor()
	cursor.execute("select end-start+1 as L from sets_for_nodes")
	Ls = map(lambda x: x[0], cursor.fetchall())
	X = range(0,10**4,20)#[0., 100., 500., 1000., 2000., 5000., 10000.]
	Y = histogram(Ls, X)[0]
	max_y_at = Y.argmax()
	msg = "Max count between {0}-{1}".format(X[max_y_at],X[max_y_at+1])
	draw_histogram(X[1:], Y, {'filename': 'test.png',\
			'X-axis': 'length', \
			'Y-axis': 'count', \
			'XA': 0., 'XE': 2000., 'XOR': 0., 'XSTP': 500., \
			'YA': 0., 'YE': max(Y), 'YOR': 0., 'YSTP': 2000., \
			'title': ['ALLFirm counts of node lengths', msg]})
	conn.close()

def histogram_of_IGR_lengths(fasta_filename):
	sizes = count_fasta_file_sizes(fasta_filename)
	X = range(0,10**4,20)#[0., 100., 500., 1000., 2000., 5000., 10000.]
	Y = histogram(sizes, X)[0]
	max_y_at = Y.argmax()
	msg = "Max count between {0}-{1}".format(X[max_y_at],X[max_y_at+1])
	draw_histogram(X[1:], Y, {'filename': 'test.png',\
			'X-axis': 'length', \
			'Y-axis': 'count', \
			'XA': 0., 'XE': 2000., 'XOR': 0., 'XSTP': 500., \
			'YA': 0., 'YE': max(Y), 'YOR': 0., 'YSTP': 2000., \
			'title': ['ALLFirm counts of IGR lengths', msg]})

def histogram_of(input, title):	
	X = range(0,10**4,20)
	Y = histogram(input, X)[0]
	max_y_at = Y.argmax()
	msg = "Max count between {0}-{1}".format(X[max_y_at],X[max_y_at+1])
	draw_histogram(X[1:], Y, {'filename': 'test.png',\
			'X-axis': 'length', \
			'Y-axis': 'count', \
			'XA': 0., 'XE': 2000., 'XOR': 0., 'XSTP': 500., \
			'YA': 0., 'YE': max(Y), 'YOR': 0., 'YSTP': 2000., \
			'title': [title, msg]})
	conn.close()

def draw_histogram(X, Y, settings):
	n = len(X)
	dislin.metafl ('png')
	dislin.setfil (settings['filename'])
	dislin.disini()
	dislin.complx()
	dislin.pagera()
	dislin.axspos (450, 1800)
	dislin.axslen (2200, 1200)
	dislin.name (settings['X-axis'], 'X')
	dislin.name (settings['Y-axis'], 'Y')
	dislin.labdig (-1, 'X')
	dislin.labdig (-1, 'Y')
	dislin.ticks (9, 'X')
	dislin.ticks (10, 'Y')
	for i,t in enumerate(settings['title']):
		dislin.titlin (t, i+1)
	ic = dislin.intrgb (0.95, 0.95, 0.95)
	dislin.axsbgd (ic)
	dislin.graf (settings['XA'], settings['XE'], settings['XOR'], settings['XSTP'],\
		     settings['YA'], settings['YE'], settings['YOR'], settings['YSTP'])
	dislin.setrgb (.7, .7, .7)
	dislin.color ('fore')
	dislin.height (50)
	dislin.title ()
	dislin.color ('blue')
	dislin.curve (X, Y, n)
	dislin.disfin ()

def count_unique_accessions(IGR_head_filename):
	acc = set()
	with open(IGR_head_filename) as f:
		for line in f:
			line = line.strip()[1:]
			acc.add(line[:line.find('/')])
	return acc

def count_fasta_file_sizes(fasta_filename):
	sizes = []
	for r in SeqIO.parse(open(fasta_filename), 'fasta'):
		insort(sizes, len(r.seq))
	return sizes	

def export_cliques_fasta(IGR_filename, pickle_filename, output_dir, prefix, output_filename, SWparam=(+8,-7,-16,-2)):
	"""
		Create fasta files for each clique as preparation of running CMFinder.
		Each clique has it's own folder under <output_dir>/<size>_<dummy_count>/.
		(Careful not to overwrite the files!)
	"""
	from miscPairwiseAlign import SWnotrace
	conn = CONN_FUNC()
	cursor = get_dict_cursor(conn)
	reader = FastaReader(IGR_filename, dun_parse_id=True)

	FETCH_SQL = "SELECT n.id,s.start,s.end \
				FROM sets_for_nodes s \
				LEFT JOIN nodes_to_index AS n \
				ON (s.nodes_ind=n.ind) WHERE i={i}"

	raw_input("Sure connection " + str(CONN_FUNC) + " is right?")
	raw_input("Did you remember to restore the DB to full set of parsed?")
	raw_input("Is the clique pickle a list of [(dummy_index, [i1,i2,i3......])]?")
	raw_input("Checked SW params? " + str(SWparam))
	raw_input("Sure the directory " + str(output_dir) + " to write to is clean?")

	with open(output_filename, 'w') as cmd_f:
#		dummy_count = 0
		with open(pickle_filename) as handle:
			QQQ = load(handle)
		for dummy_count,Q in QQQ:
			print >> sys.stderr, "clique is", Q
			if len(Q) < CLIQUE_MIN_SIZE:
				continue
#			dummy_count += 1
			tmp_name = "{prefix}_{dummy}_size{size}".format(prefix=prefix,size=len(Q),dummy=dummy_count)
			tmp_dir  = os.path.join(output_dir, tmp_name)
			cmd_f.write("cmfinder.pl -hmm -combine {0}/{1}.fna\n".format(tmp_dir, tmp_name))
			os.mkdir(tmp_dir)
			with open(os.path.join(tmp_dir, tmp_name+'.fna'), 'w') as fasta_f:
				cursor.execute( FETCH_SQL.format(i=Q[0]) )
				r = cursor.fetchone()
				acc,strand,start,end = parsed_accID(r['id'],False,r['start'],r['end'])
				seq1 = reader[r['id']].seq[r['start']-1:r['end']]
				cursor.execute( FETCH_SQL.format(i=Q[1]) )
				r = cursor.fetchone()

#				seq2 = reader[r['id']].seq[r['start']-1:r['end']]
#				r_seq1 = seq1.reverse_complement()
#				print >> sys.stderr, "Running SWnotrace #1....seq len {0},{1}".format(len(seq1),len(seq2))
#				conf1 = SWnotrace(*SWparam,s1=seq1,s2=seq2)['score']
#				print >> sys.stderr, "Running SWnotrace #2...."
#				conf2 = SWnotrace(*SWparam,s1=r_seq1,s2=seq2)['score']
#				if conf1 > conf2:
#					print >> sys.stderr, "conf1"
#					fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,seq1))
#					acc,strand,start,end = parsed_accID(r['id'],False,r['start'],r['end'])
#					fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,seq2))
#				else:
#					print >> sys.stderr, "conf2"
#					fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,r_seq1))
#					acc,strand,start,end = parsed_accID(r['id'],False,r['start'],r['end'])
#					fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,seq2))
#					seq1 = r_seq1
	
				# instead of using the time-consuming SWnotrace, re-grab the clique from parsed
				tmp = ",".join( map(str, Q) )
				cursor.execute("select i1,i2,opposite_strand "
						"from parsed "\
						"where i1 in ({0}) and i2 in ({0})".format( tmp ))
				G = Graph()
				for r in cursor.fetchall():
					if r['opposite_strand'] == 0:
						# means same-strand
						G.add_edge(r['i1'],r['i2'],+1)
					else:
						G.add_edge(r['i1'],r['i2'],-1)

				if strand == +1:
					fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,seq1))
				else:
					fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,end,start,seq1))

#				for i in Q[2:]:
				for i in Q[1:]:
					print >> sys.stderr, "fetching for i={0}....".format(i)
					cursor.execute( FETCH_SQL.format(i=i) )
					r = cursor.fetchone()
					acc,strand,start,end = parsed_accID(r['id'],False,r['start'],r['end'])
					seq = reader[r['id']].seq[r['start']-1:r['end']]
					r_seq = seq.reverse_complement()

					if G.has_edge( Q[0], i):
						direction = G.get_edge( Q[0], i )
						print >> sys.stderr, "direct edge {0} -- {1}, direction {2}".format(\
								Q[0], i, direction)
					else:
						# there's no direct edge from Q[0] -- i
						# so there must be some x, Q[0] -- x -- i
						for x in G.neighbors_iter( i ):
							if G.has_edge( Q[0], x ):
								direction = G.get_edge( Q[0], x) * G.get_edge( x, i)
								print >> sys.stderr, "indirect edge {0} -- {1} -- {2}, direction {3}".format(Q[0], x, i, direction)
								break

#					if SWnotrace(*SWparam,s1=seq1,s2=seq)['score'] > SWnotrace(*SWparam,s1=seq1,s2=r_seq)['score']:
					if direction == +1:
						if strand == +1:
							fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,seq))
						else:
							fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,end,start,seq))
					else:
						if strand == +1:
							fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,end,start,r_seq))
						else:
							fasta_f.write(">{0}/{1}-{2}\n{3}\n".format(acc,start,end,r_seq))
#							print >> sys.stderr, "check out"
#							print >> sys.stderr, ">{0}/{1}-{2}\n{3}\n".format(acc,start,end,r_seq)
#							raw_input()

if __name__=="__main__":
	#export_cliques_fasta(IGR_filename, pickle_filename, output_dir, prefix, output_filename, SWparam=(+8,-7,-16,-2))
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-i", "--igr", help="IGR filename")
	parser.add_option("-c", "--clique", help="clique pickle filename")
	parser.add_option("-d", "--output-dir", dest="output_dir", help="output directory")
	parser.add_option("-p", "--prefix", help="output motif directory prefix")
	parser.add_option("-o", "--output-filename", dest="output_filename", \
			help="output command filename", default="cmd")
	parser.add_option("-w", "--SWparam", default=(+8,-7,-16,-2), \
			help="BLAST parameters (match,mismatch,gap open,gap extend) used")

	(options, args) = parser.parse_args()
	export_cliques_fasta(options.igr, \
			options.clique, \
			options.output_dir, \
			options.prefix, \
			options.output_filename, \
			options.SWparam)
