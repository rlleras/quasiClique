from cPickle import *
import os,re,sys,fnmatch
a = load(open('ALLFirm_RefSeq25_m30s0_M8N7Q16R2W3E2_cut35.quasi80_ALL.total_motif_ranks.pickle'))
aa = filter(lambda x: x[1] is not None, a)
for x in aa:
	try:
		i = x[3].index('.fna')
		dir = x[3][:i]
	except:
		print >> sys.stderr, x[3]
		i = x[3].index('_split')
		dir = x[3][:i]
	for motif in fnmatch.filter(os.listdir('tmp/motif_dir/' + dir), '*.motif.*.dup_rmed'):
		print("""condor_run "cmbuild /4/lachesis/MOTIF_storage/ALLFirm_RefSeq25_m30s0_M8N7Q16R2W3E2/{dir}/{file}.cmbuilded /4/lachesis/MOTIF_storage/ALLFirm_RefSeq25_m30s0_M8N7Q16R2W3E2/{dir}/{file}; cmsearch --fil-T-hmm 10 /4/lachesis/MOTIF_storage/ALLFirm_RefSeq25_m30s0_M8N7Q16R2W3E2/{dir}/{file}.cmbuilded /4/lachesis/IGR_storage/ALL_knownRibo_20071102_and_randnonribo.fna > /4/lachesis/MOTIF_storage/ALLFirm_RefSeq25_m30s0_M8N7Q16R2W3E2/{dir}/{file}.hmmT10.cmscanned" &""".format(dir=dir, file=motif))
