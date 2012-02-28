import os,re,sys
from collections import defaultdict

input_real_prefix = 'output/input_fnas/ALLActinobacteria_20081002_plus20071102_curated_withDup/ALLActinobacteria_20081002_plus20071102_curated_withDup'
output_blast_prefix = 'output/output_blast_etc/ALLActinobacteria_20081002_plus20071102_curated_withDup/ALLActinobacteria_20081002_plus20071102_curated_withDup'

blast_cmd_patterns = [\
			'M5N4Q10R6W3E10.plus_self10_modishuffled',\
			'M5N4Q10R6W3E10.plus_self60_modishuffled',\
			'M5N4Q10R6W3E10.plus_prop60_modishuffled'\
		     ]
titles = ['self10','self60','prop60']
colors = [3,-1,4,5,9,2]

if 1:
	family = sys.argv[1]

	plot_cmd = 'plot '
	label_cmds = []
	FDR_ticks = [.05,.1,.15,.25,.50]
	label_spacing = 0.015
	dummy = 1
	for p,t in zip(blast_cmd_patterns,titles):
		per_FDR = defaultdict(lambda: {'score':None,'TDR':0.})

		data = "{0}_{1}.fna.{2}.WUblast".format(output_blast_prefix,family,p)
		# get timing file
		total_runtime = int(round(float(os.popen("tail -n 1 {0}.timing".format(data)).readline().split()[1])))
#		os.system("python evaluate_blast_result.py {1}_{0}.fna {2}_{0}.fna.{3}.WUblast {0} > {4}.ROC".format(family,input_real_prefix,output_blast_prefix,p,data))
		c = colors.pop(0)
		colors.append(c)
		plot_cmd += " \"{0}.ROC\" t \"{1}, {2}s\" with lines lw 3 lt {3},".format(data,t,total_runtime,c)
		with open(data+'.ROC') as f:
			for line in f:
				FDR,TDR,score_cutoff = map(float, line.strip().split())
				score_cutoff = int(round(score_cutoff))
				for x in FDR_ticks:
					if FDR <= x and TDR > per_FDR[x]['TDR']:
						per_FDR[x] = {'score':score_cutoff,'TDR':TDR,'FDR':FDR}
		for k,v in per_FDR.iteritems():		
			if v['score'] is not None:
				label_cmds.append("set label \"{0}\" at first {1}, first {2} tc lt {3}".format(v['score'],v['FDR']+label_spacing,v['TDR']-0.01,c))
		dummy += 1
        import Gnuplot
        p = Gnuplot.Gnuplot()
        p('set terminal png')
#       p('set xra [0:1]')
#       p('set yra [0:1]')
        p("set out \"test{0}.png\"".format(family))
	p("set title \"Actino {0} +5/-4,W3\"".format(family))
#	p('set title "Actino family vs family+self-dishuffled-10 +5/-4,10/6,W3"')
#	p('set xlabel "(query_length + target_length)/2"')
	p('set xlabel \"(1 - specificity)\"; set ylabel \"sensitivity\"')
#	p('set ylabel "maximal observed hit score"')
	print 'label_cmds:',label_cmds
	for cmd in label_cmds:
		p(cmd)
	p(plot_cmd[:-1])
#        p("plot \"{0}\" t \"{2}\" with lines lw 5, \"{1}\" t \"{3}\" with lines lw 1".format(data1,data2,t1,t2))
#	p('plot "ALLActinobacteria_20081002_plus20071102_curated_withDup_ALL.fna.M5N4Q10R6W3E10.plus_self10_modishuffled.WUblast.avglen_mine_scatterplot" with points lt -1 ')
        p('set out')

