import os,re,sys,fnmatch

VALID_LEFT = ['(', '<', '[']
VALID_RIGHT = [')', '>', ']']

def manual_handle():
	for file in fnmatch.filter(os.listdir('.'), '*.1.motif.*'):
		os.system("rm_dup_instance.pl {0} {0}.dup_rmed".format(file))
		os.system("rank_cmfinder.pl -w -rank {0}.dup_rmed {0}.dup_rmed.summary".format(file))
	combine_rank_summary()

def combine_rank_summary(suffix='.summary'):
	"""
	Combine all the .summary files (or .<suffix> files) in current directory
	"""
	summaries = []
	for x in fnmatch.filter(os.listdir('.'), '*'+suffix):
		print >> sys.stderr, "match ", x
		f = open(x)
		header = f.readline() # first line just the header
		for line in f:
			line = line.strip()
			if len(line) == 0:
				break
			summaries.append( (float(line.split(',')[-1]), line) )
		f.close()
	summaries.sort()
	summaries.reverse()
	print(header + "\n".join(map(lambda x: x[1], summaries)))

def furnish_motif(cons, rf):
	"""
	Given the consensus SS_cons and RF from CMFinder's motif file,
	furnish it so that cons is just ( ) .
	and rf is all upper case with gap positions removed
	"""
	assert len(cons) == len(rf)
	new_cons = ''
	new_rf = ''
	for i in xrange(len(cons)):
		if rf[i]=='.':
			continue
		if cons[i] in VALID_LEFT:
			new_cons += '('
		elif cons[i] in VALID_RIGHT:
			new_cons += ')'
		else:
			new_cons += '.'
		new_rf += rf[i].upper()
	return new_cons, new_rf


def read_motif(motif_filename):
	SS_cons = ''
	RF = ''
	with open(motif_filename) as f:
		for line in f:
			if line.startswith('#=GC SS_cons'):
				SS_cons += line.strip().split()[2]
			elif line.startswith('#=GC RF'):
				RF += line.strip().split()[2]
	return SS_cons,RF

def combine_and_rm_dup(motif_dir, only_rm_dup=False):
	"""
		Uses CMFinder_03 to combine motifs and remove duplicates.
		Doesn't really execute the commands, just prints out
		the commands so it's safer.
	"""
	old_dir = os.popen("pwd").read().strip()
	for dir in os.listdir(motif_dir):
		if not os.path.isdir(motif_dir+dir):
			continue
		print("cd {0}".format(motif_dir+dir))
#		motif_files = set(fnmatch.filter(os.listdir(motif_dir+dir),'*motif*'))
#		motif_files = motif_files.difference(set( fnmatch.filter(motif_files,'*cm*motif*') ))
#		motif_files = " ".join(motif_files)
#		if not only_rm_dup:
#			print("comb_motif.pl {0}.fna {1}".format(dir, motif_files))
		print("rm *.r *.dup_rmed *.summary")
		print("rm_dup.pl {0}.fna \"{0}.fna*.1.motif.*\"".format(dir))
		if not only_rm_dup:
			print("comb_motif.pl {0}.fna \"{0}.fna*.1.motif.*\"".format(dir))
		print("rm_dup.pl {0}.fna \"{0}.fna*.1.motif.*\"".format(dir))
		print("cd " + old_dir)

def cleaning_orphan_cm(motif_dir):
	"""
		After merging and removing redundant motifs,
		the cm partner of the removed motifs will be 
		left there as orphans, so we use this to clean
		them out since rm_dup.pl doesn't remove corressponding
		cm files.

		NOTE: keeping the orphan cm files around isn't that
		harmful, but accidentally killing useful cm files
		ARE! so don't use this unless really really sure.

		Output: a list of `rm xxxx` commands for deleting cms.
	"""
	if motif_dir is None:
		for cm_file in fnmatch.filter(os.listdir('.'),'*.cm.*'):
			if cm_file.find('cmsearch') > 0: continue # that was a cmsearch scan file, not cm file
			motif_file = cm_file.replace('.cm.','.motif.')
			if not os.path.exists(motif_file):
				print("rm " + cm_file)
		return

	for dir in os.listdir(motif_dir):
		if not os.path.isdir(dir): continue
		for cm_file in fnmatch.filter(os.listdir(motif_dir+dir),'*.cm.*'):
			if cm_file.find('cmsearch') > 0: continue # that was a cmsearch scan file, not cm file
			motif_file = cm_file.replace('.cm.','.motif.')
			if not os.path.exists(motif_dir+dir+'/'+motif_file):
				print "rm %s" % (motif_dir+dir+'/'+cm_file)	
		#print "rm %s*temp*" % (motif_dir+dir+'/')


if __name__=="__main__":
	combine_and_rm_dup('MOTIF_DIR/ALLFirm_M5N4Q8R6_motif_dir/')
	#combine_rank_summary()
	#cleaning_orphan_cm(None)
	#manual_handle()
