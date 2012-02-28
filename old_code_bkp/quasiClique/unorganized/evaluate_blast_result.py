import os,re,sys,fnmatch
from collections import defaultdict
from miscParses import parseWUBLASTline,parsed_accID,calc_overlap
from functools import reduce
from networkx import *
from networkx.readwrite import *
from Bio import SeqIO
from interval import *
from pClique import *
from miscBio import NewickReader

def parse_id(id, no_padding=True):
	"""
		If it's not padded, then it should be in in form <family>_<db_id>
		Otherwise, <family>_<db_id>_<acc>/<start>-<end>
	"""
	rex1 = re.compile('(\S+)_(\d+)')
	rex2 = re.compile('(\S+)_(\d+)_(\S+)/(\d+)-(\d+)')

	if no_padding:
		m = rex1.match(id)
		return m.group(1),m.group(2)
	else:
		m = rex2.match(id)
		return m.group(1),m.group(2),m.group(3),int(m.group(4)),int(m.group(5))

def blast_to_graph(real_fna_input,blast_output,score_cutoff,multi_edge,program='WU'):
	"""
		Make a node for each seq in <real_fna_input> (so we'll know loner nodes)
		Read <blast_output>, for each hit >= score_cutoff,
 		  make an edge (n1,n2) where n1,n2 are query & target seq IDs,

		If multi_edge is False, then
		  the edge weight is the (maximum observed) hit score
		Otherwise,
 		  the edge(s) are (hit_local_start,hit_local_end,hit_score)
		
		Stores the networkx XGraph in a pickle <blast_output>_cut<score_cutoff>.gpickle
		Then returns the pickle filename
	"""
	if program!='WU': raise Exception, 'temporarily does not support non-WUBLAST output!(TODO)'	
	
	X = XGraph()
	X.ban_selfloops() # this is for ignoring self-hits
	if multi_edge:
		X.allow_multiedges()
	else:
		X.ban_multiedges()

	with open(real_fna_input) as f:
		for r in SeqIO.parse(f,'fasta'):
			X.add_node(r.id)

	with open(blast_output) as f:
		for line in f:
			hit_dict = parseWUBLASTline(line)
			if hit_dict['sprime'] < score_cutoff: continue
			if multi_edge:
				if hit_dict['id1'] < hit_dict['id2']:
					X.add_edge(hit_dict['id1'],hit_dict['id2'],((hit_dict['start1'],hit_dict['end1']),(hit_dict['start2'],hit_dict['end2']),hit_dict['sprime']))
				else:
					X.add_edge(hit_dict['id1'],hit_dict['id2'],((hit_dict['start2'],hit_dict['end2']),(hit_dict['start1'],hit_dict['end1']),hit_dict['sprime']))
			elif (not X.has_edge(hit_dict['id1'],hit_dict['id2'])) or X.get_edge(hit_dict['id1'],hit_dict['id2']) < hit_dict['sprime']:
				X.add_edge(hit_dict['id1'],hit_dict['id2'],hit_dict['sprime'])

		
	gpickle_filename = blast_output+'_cut'+str(score_cutoff)+'.gpickle'
	write_gpickle(X,gpickle_filename)
	return gpickle_filename

def evaluate_blast_graph(graph_or_filename,ignore_prefix=['shuffled_','random_']):
	"""
		Reads the blast-processed XGraph (or the filename of the pickled XGraph)
		The REAL ncRNA seq IDs should be in format <family>_<DB id>
		The FAKE seq IDs should be in format <ignore_prefix>xxxxx.....
	
		Returns (result,sens_by_family,spec_by_family)
		where result is a dict with seqID --> (# of TP, out degree)
		sens_by_family is dict with family --> list of sens(ratio) for all nodes in this family
		spec_by_family is dict with family --> list of spec(ratio) for all nodes in this family
	"""
	if type(graph_or_filename) is XGraph:
		X = graph_or_filename
	else:
		X = read_gpickle(graph_or_filename)
	
	result = {}
	spec_by_family = defaultdict(lambda: [])
	sens_by_family = defaultdict(lambda: [])
	total_size_by_family = defaultdict(lambda: 0)
	for n in X.nodes_iter():
		i = n.rfind('_')
		family,db_id = n[:i],n[(i+1):]
		if any(map(lambda x: family.startswith(x), ignore_prefix)): continue
		total_size_by_family[family] += 1
		neighbors = X.neighbors(n)
		if len(neighbors) == 0:
			result[n] = 0
			spec_by_family[family].append(0)
			sens_by_family[family].append(0)
			continue
		good_count = reduce(lambda acc,x: acc+(family==x[:x.rfind('_')]), neighbors, 0)
		result[n] = (good_count,len(neighbors))
		spec_by_family[family].append(good_count*1./len(neighbors))
		sens_by_family[family].append(good_count) # is raw count!!

	for k in sens_by_family:
		sens_by_family[k] = map(lambda x: x*1./total_size_by_family[k], sens_by_family[k])
	return (result,sens_by_family,spec_by_family)

def evaluate_blast_graph_count_nucleotides(graph_or_filename,hit_ratio=None,ignore_prefix=['shuffled','random']):
	"""
		Similar to evaluate_blast_graph, except that the real ncRNAs (in the DB, not query)
		  are embedded with flanking regions, and the IDs should be in 
		  format <family>_<DB id>_<acc>/<embedded_start>-<embedded_end>

		If <hit_ratio> is None, then for each node N, 
		  sensitivity = (# of real ncRNA-neighbor nts) / (# of real ncRNA nts)
		  specificity = (# of real ncRNA-neighbor nts) / (# of neighbor nts)
		
		If <hit_ratio> is defined, ex: 0.8, then for each node N,
     		  a neighbor node M is a hit if the # of hit on M is >= <hit_ratio>*<M's ncRNA len>

		NOTE: for this kind of blast output, the INPUT should be seq IDs like <family>_<db_id>
		      which means they are real ncRNAs with NO padding
		      and the DB can either be random/shuffled seqIDs
 		      or <family>_<db_id>_<acc>/<embedded_start>-<embedded_end>
	"""
	rex = re.compile('(\S+)_(\d+)_(\S+)')
	rex_real = re.compile('(\S+)_(\d+)')
	from miscncRNA import get_ncRNA_info

	if hit_ratio is not None:
		hit_ratio = float(hit_ratio)

        if type(graph_or_filename) is XGraph:
                X = graph_or_filename
        else:
                X = read_gpickle(graph_or_filename)

	total_nt_by_family = defaultdict(lambda: 0)
        spec_by_family = defaultdict(lambda: [])
        sens_by_family = defaultdict(lambda: [])
	for n in X.nodes_iter():
		if any(map(lambda x: n.startswith(x), ignore_prefix)): continue
		if n.count('_') > 2: continue
		m = rex_real.match(n)
		if m is None: continue
		family,query_db_id = m.group(1),m.group(2)
		tmp_true = defaultdict(lambda: IntervalSet())
		tmp_false = defaultdict(lambda: IntervalSet())

		# the query nodes must be <family>_<db_id> (i.e. no padding)
		info = get_ncRNA_info(query_db_id)
		print >> sys.stderr, n
		total_nt_by_family[family] += info['end']-info['start']+1
		if X.degree(n) == 0: # has 0 neighbors
			sens_by_family[family].append(0)
			spec_by_family[family].append(0)
			continue		
		# e is in format (local_start,local_end,score)
		for (myself,neighbor,e) in X.edges_iter(n):
			if any(map(lambda x: neighbor.startswith(x), ignore_prefix)):
				# not a real ncRNA
				tmp_false[neighbor].add(Interval(e[0],e[1]))
			else:
				m = rex.match(neighbor)
				duncare,db_id,blob = m.group(1),m.group(2),m.group(3)
				if db_id == query_db_id: continue # it's a self vs self-embedded hit, ignore
				(acc,duncare),hit_start,hit_end,hit_strand = parsed_accID(blob,True,e[0],e[1])
				tmp_true[db_id].add(Interval(hit_start,hit_end))

		tp,fp = (0,0)
		for db_id,regions in tmp_true.iteritems():
			info = get_ncRNA_info(db_id)
			for x in regions:
				c = calc_overlap(info['start'],info['end'],x.lower_bound,x.upper_bound)
				if hit_ratio is None:
					tp += c
					fp += (x.upper_bound-x.lower_bound+1) - c
					
				elif c >= hit_ratio*(info['end']-info['start']+1):
					tp += 1
				else:
					fp += 1
		for some_id,regions in tmp_false.iteritems():
			for x in regions: fp += x.upper_bound-x.lower_bound+1

		print >> sys.stderr, tp,fp
		if tp+fp == 0:
			sens_by_family[family].append(0)
			spec_by_family[family].append(0)
		else:
			sens_by_family[family].append(tp) # NOTE: it's raw count!!!
			spec_by_family[family].append(tp*1./(tp+fp))
		#raw_input('...')
	for k in sens_by_family:
		if hit_ratio is None:
			sens_by_family[k] = map(lambda x: x*1./total_nt_by_family[k], sens_by_family[k])
		else:
			sens_by_family[k] = map(lambda x: x*1./len(total_nt_by_family[k]), sens_by_family[k])
	return (None,sens_by_family,spec_by_family)

def evaluate_blast_graph_by_phylo(graph_or_filename,phylo_filename,ignore_prefix=['shuffled_','random_']):
	"""
		Reads the blast-processed XGraph (or the filename of the pickled XGraph)
		The REAL ncRNA seq IDs should be in format <family>_<DB id>
		The FAKE seq IDs should be in format <ignore_prefix>xxxxx.....
	
		Returns (result,sens_by_family,spec_by_family)
		where result is an obsolete junk (so just make it None) for now
		phylo_sum_by_family is dict with family --> list of sums of neighbor phylo distances 
							    for all nodes in this family
		spec_by_family is dict with family --> list of spec(ratio) for all nodes in this family
	"""
	if type(graph_or_filename) is XGraph:
		X = graph_or_filename
	else:
		X = read_gpickle(graph_or_filename)

	p = NewickReader(phylo_filename)	
	
	phylo_sum_by_family = defaultdict(lambda: [])
	spec_by_family = defaultdict(lambda: [])
	total_size_by_family = defaultdict(lambda: 0)
	for n in X.nodes_iter():
		i = n.rfind('_')
		family,db_id = n[:i],n[(i+1):]
		if any(map(lambda x: family.startswith(x), ignore_prefix)): continue
		total_size_by_family[family] += 1
		if X.degree(n) == 0:
			spec_by_family[family].append(0)
			phylo_sum_by_family[family].append(0)
			continue
		phylo_sum,tp = 0,0
		for m in X.neighbors_iter(n):
			if any(map(lambda x: family.startswith(x), ignore_prefix)): continue
			m1 = parse_id(m,True)			
			if m1[0] == family:
				try:
					phylo_sum += p.distance(db_id,m1[1])
				except:
					pass
				tp += 1	
		phylo_sum_by_family[family].append(phylo_sum)		
		spec_by_family[family].append(tp*1./X.degree(n))

	return (None,phylo_sum_by_family,spec_by_family)

def plotROC_for_family(real_fna_input,blast_output,family,count_nucleotides,phylo_filename=None):
	x = []
	y = []
	s = []
	start_score_cutoff,stop_score_cutoff,step_score_cutoff = (10,1000,2)

	X = read_gpickle(blast_to_graph(real_fna_input,blast_output,start_score_cutoff,count_nucleotides))

	last_num_of_edges = -1

	for score_cutoff in xrange(start_score_cutoff, stop_score_cutoff, step_score_cutoff):
		XX = X.copy()
		if count_nucleotides:
			XX.delete_edges_from(filter(lambda e: float(e[2][2]) < score_cutoff, XX.edges_iter()))
		else:
			XX.delete_edges_from(filter(lambda e: e[2] < score_cutoff, XX.edges_iter()))

		if XX.number_of_edges() == 0: break
		if XX.number_of_edges()==last_num_of_edges: continue # no need to run it since results are same
		last_num_of_edges = XX.number_of_edges()
		print >> sys.stderr, "{0} edges, {1} nodes with cutoff {2}....".format(XX.number_of_edges(),XX.number_of_nodes(),score_cutoff)

		# TEST, DELETE LATER
		#(m_cliques, id_map) = maximal_cliques(XX)
		#print id_map
		#print m_cliques
		#sys.exit(-1)
		
		if count_nucleotides:
			(result,sens_by_family,spec_by_family) = evaluate_blast_graph_count_nucleotides(XX)#evaluate_blast_graph(XX)
		else:
			(result,sens_by_family,spec_by_family) = evaluate_blast_graph_by_phylo(XX,phylo_filename)

		print >> sys.stderr, "sanity check {0} sens_by_family[{1}], {2} sens_by_family[{3}]".format(len(sens_by_family[family]),family,len(spec_by_family[family]),family)
#		print result
		#print >> sys.stderr, spec_by_family[family]
		if spec_by_family[family].count(0)==len(spec_by_family[family]): break
		# calculate the average sensitivity and specificity
		sens = sum(sens_by_family[family])*1./max(1,len(sens_by_family[family])-sens_by_family[family].count(0))
		spec = sum(spec_by_family[family])*1./max(1,len(spec_by_family[family])-spec_by_family[family].count(0))
		print >> sys.stderr, "avg sens: {0}, avg spec: {1}".format(sens,spec)
		x.append(1-spec)
		y.append(sens)
		s.append(score_cutoff)
	for i in xrange(len(x)):
		print x[i],y[i],s[i]
	return x,y,s	

def scatterplot_length_to_maxscore(dirname,filename_pattern):
	"""
		Read through ALL files matching <filename_pattern> description
		under directory <dirname>

		Assumes the files are in WU-BLAST format		
		Then plots a scatterplot where 
		  x-axis is the length(query+target seq)
		  y-axis is the max score observed for SOME pairs of seq where sum of length is x
	"""
	maxscore_of_length = defaultdict(lambda: 10)
	for filename in fnmatch.filter(os.listdir(dirname),filename_pattern):
		print >> sys.stderr, filename
		with open(dirname+'/'+filename) as f:
			for line in f:
				hit_dict = parseWUBLASTline(line)				
				l1 = abs(hit_dict['end2']-hit_dict['start2'])+1
				l2 = abs(hit_dict['end1']-hit_dict['start1'])+1
				l = (l1+l2)/2
				maxscore_of_length[l] = min(maxscore_of_length[l], hit_dict['e'])
	for l,s in maxscore_of_length.iteritems():
		print l,s

def scatterplot_phylodist_to_score(blast_dirname,filename_rex,phylo_dirname,NONMEMBER_DIST=10.):
	"""
		Similar to scatterplot_length_to_maxscore, except that we're plotting
		the branch distance to maximum observed score. And we do it PER family
		blast file since branch distances make sense only within a family....

		Also:
		  We assume that the blast files not only match filename_rex, but also
		  ends in '.WUblast'

		  We assume that under <phylo_dirname>, there will only be 1 .outtree file
		  that fits the filename_rex-extracted family name

		  We also assume that this blast output is a "real" vs "padded_real + negative",
		  so id1 is always not padded, and id2 is potentially padded
	"""
	for filename in os.listdir(blast_dirname):
		if os.path.exists(blast_dirname+'/'+filename+'.phylodist_to_score'):
			continue

		maxscore_of_dist = defaultdict(lambda: float("-inf"))
		m = filename_rex.match(filename)
		if m is None or not filename.endswith('WUblast'):
			continue
		family = m.group(1)
		phylo_filename = fnmatch.filter(os.listdir(phylo_dirname),"*{0}*.outtree".format(family))[0]
		p = NewickReader(phylo_dirname+'/'+phylo_filename)
		print >> sys.stderr, filename
		with open(blast_dirname+'/'+filename,'r') as f:
			for line in f:
				d = parseWUBLASTline(line)
				id1,id2 = d['id1'],d['id2']
				if id1[:id1.find('_')]==id2[:id2.find('_')] and id1.startswith(family):
					try:
						dist = p.distance(parse_id(id1,True)[1], parse_id(id2,False)[1])
						maxscore_of_dist[dist] = max(maxscore_of_dist[dist], d['sprime'])
					except:
						print >> sys.stderr, "failed phylo dist on ",id1,id2
				else: # a non-hit
					maxscore_of_dist[NONMEMBER_DIST] = max(maxscore_of_dist[NONMEMBER_DIST], d['sprime'])

		dists = maxscore_of_dist.keys()
		dists.sort()
		with open(blast_dirname+'/'+filename+'.phylodist_to_score','w') as f:
			for k in dists:
				f.write("{0}\t{1}\n".format(k, maxscore_of_dist[k]))

class Drawer:
	def __init__(self):
		self.title = None
		self.xlabel = None
		self.ylabel = None
		self.plotstyle = None
		self.output_png_name = None
		self.plot_cmd = ''

	def add_plot(self,data,legend,c,pts):
		self.plot_cmd += " \"{0}\" t \"{1}\" with {4} lw 1 lt {2} pt {3},".format(data,legend,c,pts,self.plotstyle)

	def clear_plot_cmd(self):
		self.plot_cmd = ''

	def draw(self):
		from Gnuplot import Gnuplot
	        p = Gnuplot()
        	p('set terminal png')
		p("set size 0.60,0.60")
		p("set key right bottom")
	        p("set out \"{0}.png\"".format(self.output_png_name))
		p("set pointsize 0.6")
		p("set title \"{0}\"".format(self.title))
		p("set xlabel \"{0}\"; set ylabel \"{1}\"".format(self.xlabel, self.ylabel))
		p("plot "+self.plot_cmd[:-1])
        	p('set out')

def test_drawer():
	d = Drawer()
	#d.title = 'test'
	d.xlabel = 'phylo-distance'
	d.ylabel = 'max hit score'
	d.plotstyle = 'points'
	#d.output_png_name = 'test'
	rex = re.compile('\S+withDup_(\S+).fna\S+')
	for filename in os.listdir('.'):
		if filename.endswith('.phylodist_to_score') and filename.find('M8N7Q16R2W3E10') > 0:
			print >> sys.stderr, filename
			m = rex.match(filename)
			family = m.group(1)
			d.clear_plot_cmd()
			d.title = "Firm {0}, +8/-7, 16/2, W3".format(family)
			d.output_png_name = 'Firm_M8N7Q16R2_'+family+'.phylodist_scatterplot'
			d.add_plot(filename,family,6,3)
			d.draw()

def test_drawer2():
	from miscncRNA import get_all_rbtypes
	colors = [1,3,4,5,6]
	d = Drawer()
	d.xlabel = '1-spec'
	d.ylabel = 'average sum-of-phylo-distance'
	d.plotstyle = 'linespoints'
	rex = re.compile('\S+withDup_\S+.fna.(\S+).plus_padseqlen\S+')
	for family in get_all_rbtypes():
		d.clear_plot_cmd()
		d.title = "Firm {0}".format(family)
		d.output_png_name = 'Firm_phylo_ROC_'+family
		for filename in fnmatch.filter(os.listdir('.'),"*{0}*.phylo_ROC".format(family)):
			print >> sys.stderr, filename
			m = rex.match(filename)
			c = colors.pop(0)
			colors.append(c)
			d.add_plot(filename,family+'-'+m.group(1),c,3)
		d.draw()	

if __name__=="__main__":
	test_drawer2()
	sys.exit(-1)

#	blast_to_graph('test.fna','testM8N7Q16R2.WUblast',30,True)
#	blast_to_graph('test.fna','testM5N4Q8R6.WUblast',30,True)
	scatterplot_phylodist_to_score('output/output_blast_etc/ALLFirmicutes_20081002_plus20071102_curated_withDup_padseqlen', re.compile('ALLFirmicutes_20081002_plus20071102_curated_withDup_(\S+).fna.M8N5Q20R5W3E10.plus_padseqlen_prop100_modishuffled.WUblast'),'output/input_stos/ALLFirmicutes_20081002_plus20071102_curated/')
#	test_drawer()
	sys.exit(-1)
#	evaluate('../input_fnas/ALLFirmicutes_20081002_plus20071102_curated.fna',\
#'ALLFirmicutes_20081002_plus20071102_curated.fna.M5N4Q10R6W1E10.WUblast',40)
#	evaluate_blast_graph_count_nucleotides(blast_to_graph(sys.argv[1],sys.argv[2],float(sys.argv[3]),True))
	plotROC_for_family(sys.argv[1],sys.argv[2],sys.argv[3],True)
#	scatterplot_length_to_maxscore(sys.argv[1],sys.argv[2])
