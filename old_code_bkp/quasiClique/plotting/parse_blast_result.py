import os,re,sys
from collections import defaultdict
from miscParses import parseWUBLASTline

def extract_best_pair_hits(blast_output,program='WU'):
	"""
	For each pair of (query,target) hits, regardless of hit region length,
	take only the highest scoring pair
	"""
	if program!='WU': raise Exception, 'temporarily does not support non-WUBLAST output!(TODO)'	
	
	res = defaultdict(lambda: 0)
	binned_scores = defaultdict(lambda: 0)
	with open(blast_output) as f:
		for line in f:
			hit_dict = parseWUBLASTline(line)
			res[(hit_dict['id1'],hit_dict['id2'])] = max(res[(hit_dict['id1'],hit_dict['id2'])],hit_dict['sprime'])
	for v in res.itervalues():
		binned_scores[int(round(v))] += 1

        # now plot it
        x = binned_scores.keys()
        x.sort()
        y = [binned_scores[k] for k in x]
        for i in xrange(len(x)):
                print x[i],y[i]
        import Gnuplot
        p = Gnuplot.Gnuplot()
        p('set terminal png')
        p("set xra [0:{0}]".format(max(x)+10))
        p("set yra [0:{0}]".format(max(y)+10))
        p("set out '{0}_unique_score_dist.png'".format(blast_output))
        D = Gnuplot.Data(x,y,with_=" lines lt -1 lw 1")
        p.plot(D)
        p('set out')
	

def plot_score_distribution(blast_output,program='WU'):
	if program=='WU':
		field_pos = 5

	binned_scores = defaultdict(lambda: 0)
	for line in os.popen("cut -f {0} {1}".format(field_pos,blast_output)):
		binned_scores[int(round(float(line.strip())))] += 1

	# now plot it
	x = binned_scores.keys()
	x.sort()
	y = [binned_scores[k] for k in x]
	for i in xrange(len(x)):
		print x[i],y[i]
	import Gnuplot
	p = Gnuplot.Gnuplot()
	p('set terminal png')
	p("set xra [0:{0}]".format(max(x)+10))
	p("set yra [0:{0}]".format(max(y)+10))
	p("set out '{0}_score_dist.png'".format(blast_output))
	D = Gnuplot.Data(x,y,with_=" lines lt -1 lw 1")
	p.plot(D)
	p('set out')

if __name__=="__main__":
#	plot_score_distribution(sys.argv[1])
	extract_best_pair_hits(sys.argv[1])
#	import matplotlib
#	matplotlib.use('Agg')
