from settings import *

DIRS = {'QC': '/home/etseng/silo/UW_Larry/output/clustered/ALLFirm_sdustM8N7Q16R2_test4_part2quasi80then60_just_ncRNA_cliques',\
		'HC': '/home/etseng/silo/UW_Larry/output/clustered/ALLFirm_M8N7Q16R2_HierarchicalClustering_score50_justncRNAs'}

def compare_best_motifs(fam, criterion):
	"""
	Read the {STORE_DIR}/<fam>.txt file respectively,
	find the entry with the best <criterion> (ex: MCC or MyScore)
	and return a dict of <key>(key of DIRS) --> best entry
	"""
	result = {}
	for k,d in DIRS.iteritems():
		best = None
		filename = os.path.join(d, STORE_DIR, fam+'.txt')
		if os.path.exists(filename):
			for obj in csv.DictReader(open(filename), delimiter='\t'):
				if best is None or float(obj[criterion]) > float(best[criterion]):
					best = obj
			result[k] = best
		else:
			result[k] = defaultdict(lambda: 'NA')
		result[k]['FAMILY'] = fam
		result[k]['METHOD'] = k
	return result

if __name__ == "__main__":
	fieldnames = ['FAMILY','METHOD','MOTIF','RANK',\
			'PSCORE','PSCORE_RAxML','TP','FP','FN','TN','MCC','MyScore']
	print("\t".join(fieldnames))
	csv_out = csv.DictWriter(sys.stdout, fieldnames, delimiter='\t')
	for fam in db_summary:
		result = compare_best_motifs(fam, 'MyScore')
		for k,v in result.iteritems():
			csv_out.writerow(v)
