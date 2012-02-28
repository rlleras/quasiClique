import os,re,sys,fnmatch
from evaluate_blast_result import plotROC_for_family

d1 = 'output/input_fnas/ALLFirmicutes_20081002_plus20071102_curated_withDup'
d2 = 'output/input_stos/ALLFirmicutes_20081002_plus20071102_curated'
d3 = 'output/output_blast_etc/ALLFirmicutes_20081002_plus20071102_curated_withDup_padseqlen'

for filename in os.listdir(d1):
	family = filename[(filename.find('withDup_')+8):filename.find('.fna')]
	if family.startswith('SAM'):
		continue
	print family
	# find the sto tree file
	l = fnmatch.filter(os.listdir(d2),"*{0}*.sto.outtree".format(family))
	if len(l) == 0:
		print >> sys.stderr, "can't find phylo tree for {0}".format(family)
		continue
	phylo_filename = d2+'/'+ l[0]
	for filename2 in fnmatch.filter(os.listdir(d3),"*{0}*prop100_modishuffled.WUblast".format(family)):
		if os.path.exists(d3+'/'+filename2+'.phylo_ROC'): continue
		try:
			x,y,s = plotROC_for_family(d1+'/'+filename,d3+'/'+filename2,family,False,phylo_filename)
		except:
			continue
		f = open(d3+'/'+filename2+'.phylo_ROC','w')
		for i in xrange(len(x)):
			f.write("{0}\t{1}\t{2}\n".format(x[i],y[i],s[i]))
		f.close()
