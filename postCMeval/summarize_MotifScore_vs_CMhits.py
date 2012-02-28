from settings import *
import numpy as np

SUFFIX = '.score_vs_CMhit.evaled'
R_SCRIPT = os.path.join(os.environ['LCODE'], 'draw_boxplot_score_vs_CMhits.R')

def main(all_filename='ALL.score_vs_CMhits.evaled', delimiter='\t'):
	"""
	Given the summarized eval file which has on each line:
	FAMILY \t FAMSIZE \t MOTIF \t RANK \t PSCORE \t TOTAL_RECALL \t TPS \t FPS

	separate each line into {STORE_DIR}/<family>.{SUFFIX}
	"""
	fieldnames = ['MOTIF', 'RANK', 'PSCORE', 'TOTAL_RECALL', 'TPS', 'FPS']
	for obj in csv.DictReader(open(all_filename), delimiter=delimiter):
		fam = obj['FAMILY']
		# append this to <fam>.{SUFFIX} in {STORE_DIR}
		out = os.path.join(STORE_DIR, fam + SUFFIX)
		print >> sys.stderr, "appending to {0}....".format(out)
		os.system("echo \"{rest}\" >> {out}".format(\
				out=out,\
				rest=delimiter.join(obj[x] for x in fieldnames)))

def draw(pscore_or_rank='pscore'):
	"""
	For every .{SUFFIX} file in {STORE_DIR}, 
	sort it by rank (field 2) or pscore (field 3) 
	then draw it using {R_SCRIPT}
	"""
	pattern = os.path.join(STORE_DIR, '*'+SUFFIX)
	for file in glob.iglob(pattern):
		print >> sys.stderr, "sort then draw for {0}....".format(file)
		fam = os.path.basename(file)
		fam = fam[:-len(SUFFIX)]
		if pscore_or_rank == 'pscore':
			output = file + '.pscore_sorted'
			os.system("sort -g -k 3 {input} > {output}".format(\
					input=file, output=output))
		elif pscore_or_rank == 'rank':
			output = file + '.rank_sorted'
			os.system("sort -g -k 2 {input} > {output}".format(\
					input=file, output=output))
		else: raise Exception, "pscore_or_rank must be either 'pscore' or 'rank'!"
		os.system("R --vanilla --fam={0} --input={1} --scorename={2} < {3}".format(\
				fam, output, pscore_or_rank, R_SCRIPT))

def calc_CMscan_separation_score(tps, fps):
	"""
	My homemade scoring function for rewarding lots of TPs
	 AND TPs with much higher CM scan scores than FPs (median)

	Currently returns sum of
	      (TP_score - FP_median) / z'
	where z' is a modified z-score that 
		1) is 1. if TP_score is above average
		2) is ceil[ | z-score | ] otherwise

	So for good TP hits (higher-than-median), it contributes
	marginally w.r.t to FP_median. For bad TP hits, the
	contribution is inversely weighted by how far it is from 
	the median, so outlier TPs won't severely injure score.
	"""
	std_tp = np.std(tps) if len(tps) > 0 else 0
	std_fp = np.std(fps) if len(fps) > 0 else 0
	avg_tp = np.mean(tps) if len(tps) > 0 else 0
	avg_fp = np.mean(fps) if len(fps) > 0 else 0
	med_tp = np.median(tps) if len(tps) > 0 else 0
	med_fp = np.median(fps) if len(fps) > 0 else 0
	_tp_sum = _fp_sum = 0
	for x in tps:
		_z = 1. if x >= avg_tp else np.ceil((avg_tp - x) / std_tp)
		_tp_sum += (x - med_fp) / _z
# NOTE: not sure if I like the _fp_sum term now, ignore
#	for x in fps:
#		_z = 1. if x >= avg_fp else np.ceil((avg_fp - x) / std_fp)
#		_fp_sum += (med_tp - x) / _z
	return _tp_sum

def calc_MCC(filename, fam, output_prefix, db_summary, score_cutoff=0.):
	"""
	Read through the evaled file, and writes out to 
	<output_prefix>.txt --- (per line) 
		motif_name, rank, pscore, tp(count), fp, fn, tn, MCC, home-made-TP/FP score
	"""
	VARNA_APPLET_NAME = "VARNAv37.jar"
	COL_PER_VARNA = 2 # number of cols for Varna-applet motif drawing
	HEIGHT_PER_ROW_VARNA = 400 # per motif drawing height (px)
	WIDTH_VARNA = 1200 # width of the html page for Varna (px)
	import math
	import miscCMF

	f_html = open(output_prefix+'.html', 'w')
	f_out = open(output_prefix+'.txt', 'w')
	f_out.write("MOTIF\tRANK\tPSCORE\tTP\tFP\tFN\tTN\tMCC\tMyScore\n")

	chunk_to_write = [] # list of (MCC, chunk_dict)

	with open(filename) as f:
		i = 0
		for line in f:
			i += 1
			raw = line.strip().split('\t')
			if len(raw) == 6:
				motif_name, rank, pscore, tp_count, tps, fps = raw
			elif len(raw) == 5:
				motif_name, rank, pscore, tp_count, tps = raw
				fps = ''
			elif len(raw) == 4:
				motif_name, rank, pscore, tp_count = raw
				tps = fps = ''
			else:
				raise ValueError, "wacky!!! {0}".format(raw)
			tps = map(float, tps.split(',')) if len(tps) > 0 else []
			fps = map(float, fps.split(',')) if len(fps) > 0 else []
			TP = len(filter(lambda x: x >= score_cutoff, tps)) # the rest are FN
			FP = len(filter(lambda x: x >= score_cutoff, fps)) # the rest are TN
			FN = db_summary[fam]['TRUE'] - TP
			TN = db_summary[fam]['CONTROLS'] - FP
			MCC = ((TP * TN) - (FP * FN))*1. / max(1, math.sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)))
			myscore = calc_CMscan_separation_score(tps, fps) / max(1., TP)
			f_out.write(str(motif_name) + '\t' + str(rank)  + '\t' + str(pscore) + '\t' + str(TP) + '\t' + str(FP) + '\t' + str(FN) + '\t' + str(TN) + '\t' + str(MCC) + '\t' + str(myscore) + '\n')

			ind = motif_name.find('.')
			if ind > 0:
				motif_filename = "motifs/{0}/{1}".format(motif_name[:ind], motif_name)
			else:
				motif_filename = "motifs/{0}/{0}".format(motif_name)
			cons, rf = miscCMF.furnish_motif(*miscCMF.read_motif(motif_filename))
			chunk_to_write.append( (MCC, {'motif_name':motif_name, 'cons':cons, 'rf': rf}) )

	# sort chunk_to_write by decreasing order of MCC
	chunk_to_write.sort(key=lambda x: x[0], reverse=True)
	N = len(chunk_to_write)

	# ---------------------- VARNA APPLET HTML WRITING ---------------------- #
	rows = N / COL_PER_VARNA + (N % COL_PER_VARNA > 0)
	f_html.write("""
	<applet  code="VARNA.class"
	codebase="."
	archive="{varna}"
	width="{width}" height="{height}">
	<param name="rows" value="{rows}" />
	<param name="columns" value="{columns}" />
	""".format(\
			varna=VARNA_APPLET_NAME,\
			width=WIDTH_VARNA,\
			height=rows*HEIGHT_PER_ROW_VARNA,\
			columns=COL_PER_VARNA,\
			rows=rows\
			))
	# remember sequence/struuctureDBN<i> has to be 1-based!
	# so must i+1 when using enumerate
	for i, (MCC, chunk_dict) in enumerate(chunk_to_write):
		f_html.write("""
	<param name="sequenceDBN{i}" value="{rf}" />
	<param name="structureDBN{i}" value="{cons}" />
	<param name="titleSize{i}" value="12" />
	<param name="title{i}" value="{motif_name}(MCC {MCC:.2f})" />
	""".format(\
			i=i+1,\
			cons=chunk_dict['cons'],\
			rf=chunk_dict['rf'],\
			motif_name=chunk_dict['motif_name'],\
			MCC=MCC\
			))
	f_html.write("</applet>\n")
	# ---------------------- VARNA APPLET HTML WRITING ---------------------- #

if __name__ == "__main__":
	# ------------------------------------------------------------------
	# call main() to parse a single multi-family summary file into
	# per-family subfiles in {STORE_DIR}
	# ------------------------------------------------------------------
	main()

	# ------------------------------------------------------------------
	# call draw() after main() to draw pscore/rank vs CM hits pngs
	# ------------------------------------------------------------------
	draw('pscore')  # usage: draw('rank|pscore')
	draw('rank')
	sys.exit(-1)
	# ------------------------------------------------------------------
	# ------------------------------------------------------------------
	p_suffix = SUFFIX + '.pscore_sorted'
	for filename in glob.iglob(STORE_DIR + '/*' + p_suffix):
		fam = os.path.basename(filename)[:-len(p_suffix)]
		print >> sys.stderr, "file: {0}, family: {1}...".format(filename, fam)
		calc_MCC(filename, \
				fam=fam,\
				output_prefix=fam,\
				db_summary=db_summary)
