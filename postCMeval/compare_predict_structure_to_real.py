from settings import *
import miscCMF

"""
 Needs: description.pickle, 
"""

def main():
	handles = {}
	os.chdir(MOTIF_DIR)
	for d in os.listdir(os.path.curdir):
		fam = dir_info[d]['family']
		if fam not in handles:
			f = open(os.path.join(os.pardir, STORE_DIR, fam+'.SS_cons.compare.txt'), 'w+')
			handles[fam] = f
			handles[fam].write(">{fam}\n{ss_cons}\n".format(fam=fam, ss_cons=db_summary[fam]['SS_CONS']))
		motif_info = {}
		for file in glob.iglob(d + "/*.cmsearched"):
			motif_name = os.path.splitext(os.path.basename(file))[0]
			print >> sys.stderr, "extracting SS_cons from {0}....".format(motif_name)
			# CMFinder doesn't have the #=GC RF line, so just make a fake XXXX... string
			# since all we want is a furnished SS_cons
			ss_cons = miscCMF.read_motif(os.path.join(d, motif_name))[0]
			ss_cons = miscCMF.furnish_motif(ss_cons, 'X'*len(ss_cons))[0]
			handles[fam].write(">{motif}\n{ss_cons}\n".format(\
					motif=motif_name,\
					ss_cons=ss_cons))
			
if __name__ == "__main__":
	main()
