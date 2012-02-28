import os,re,sys
from collections import defaultdict

input_real_prefix = 'output/input_fnas/ALLFirmicutes_20081002_plus20071102_curated_withDup/ALLFirmicutes_20081002_plus20071102_curated_withDup'
output_blast_prefix = 'output/output_blast_etc/ALLFirmicutes_20081002_plus20071102_curated_withDup_padseqlen/ALLFirmicutes_20081002_plus20071102_curated_withDup'

blast_cmd_patterns = [\
#			('M5N4Q10R6W3E10.plus_self10_modishuffled',1,-1),\
#			('M5N4Q10R6W3E10.plus_self60_modishuffled',9,-1),\
#			('M5N4Q10R6W3E10.plus_prop10_modishuffled',3,6),\
#			('M5N4Q10R6W3E10.plus_prop60_modishuffled',7,6),\
#			('M5N4Q10R6W3E10.plus_prop60_random',4,3),\
#			('M5N4Q10R6W3E10.plus_random',5,3),\

#			('M5N4Q10R6W3E10.plus_padseqlen_prop60_random',1,-1),\
#			('M5N4Q10R6W11E10.plus_padseqlen_prop60_random',7,-1),\
#			('M5N4Q20R10W3E10.plus_padseqlen_prop60_random',7,3),\
			('M5N4Q8R6W3E10.plus_padseqlen_prop100_modishuffled',7,6),\
			('M8N7Q16R2W3E10.plus_padseqlen_prop100_modishuffled',5,4),\
			('M8N5Q20R5W3E10.plus_padseqlen_prop100_modishuffled',3,5),\
]
#titles = ['self-shuffled 10X','self-shuffled 60X','all-shuffled 10X','all-shuffled 60X']
#titles = ['self-shuffled 60X','all-shuffled 60X','random-nonIGR 60X']
#titles = ['self60(W3)','self60(W11)','random60(W3)','random60(W11)']
titles = ['+5/-4,8/6','+8/-7,16/2','+8/-5,20/5']
#titles = ['all-shuffled 60X','random-nonRNA 60X','random >> 1000X']

def draw_comparison_plot(family,output_png_prefix='test_'):
#	points = [1,3,5,7,6,9]
	plot_cmd = 'plot '
	label_cmds = []
	FDR_ticks = [.05,.1,.15,.25,.50]
	label_spacing = 0.015
	dummy = 1
	for (p,pts,c),t in zip(blast_cmd_patterns,titles):
		per_FDR = defaultdict(lambda: {'score':None,'TDR':0.})

		data = "{0}_{1}.fna.{2}.WUblast".format(output_blast_prefix,family,p)
		# get timing file
		total_runtime = int(round(float(os.popen("tail -n 1 {0}.timing".format(data)).readline().split()[1])))
		if not os.path.exists(data+'.ROC'):
			os.system("python evaluate_blast_result.py {1}_{0}.fna {2}_{0}.fna.{3}.WUblast {0} > {4}.ROC".format(family,input_real_prefix,output_blast_prefix,p,data))
		plot_cmd += " \"{0}.ROC\" t \"{1}, {2}s\" with linespoints lw 1 lt {3} pt {4},".format(data,t,total_runtime,c,pts)
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
	p("set size 0.60,0.60")
	p("set key right bottom")
        p("set out \"{0}_{1}.png\"".format(output_png_prefix,family))
	p("set pointsize 0.6")
	p("set title \"Firm {0} DB modishuffled 100X, W3, diff scoring\"".format(family))
#	p('set title "Firm family vs family+self-dishuffled-10 +5/-4,10/6,W3"')
#	p('set xlabel "(query_length + target_length)/2"')
	p('set xlabel \"(1 - specificity)\"; set ylabel \"sensitivity\"')
#	p('set ylabel "maximal observed hit score"')
	print 'label_cmds:',label_cmds
	for cmd in label_cmds:
		p(cmd)
	p(plot_cmd[:-1])
#        p("plot \"{0}\" t \"{2}\" with lines lw 5, \"{1}\" t \"{3}\" with lines lw 1".format(data1,data2,t1,t2))
#	p('plot "ALLFirmicutes_20081002_plus20071102_curated_withDup_ALL.fna.M5N4Q10R6W3E10.plus_self10_modishuffled.WUblast.avglen_mine_scatterplot" with points lt -1 ')
        p('set out')

def iterate_over_families(phylum,output_png_prefix):
	from miscMySQL import get_families

	for family in get_families(phylum):
		if family=='tbox': continue #TODO: delete later
		draw_comparison_plot(family,output_png_prefix)
				
if __name__=="__main__":
	iterate_over_families('Firmicutes','Firm_diffParam2')
#	draw_comparison_plot('SRP_bact','Firmg_diffGap')
