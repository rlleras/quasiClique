import os,re,sys,fnmatch
from networkx import *
from networkx.readwrite import *
from Bio import SeqIO
from miscncRNA import get_ncRNA_info
from miscParses import gen_accID

PID_GRAPH_DIR = '/homes/gws/lachesis/Larry_Riboswitch/Latest_ribo_20081002_plus20071102_curated/pid_graphs/'

def avg_pid_by_family(fasta_filename,family):
	"""
		Find the .sto.pid_graph.pickle file for <family> in PID_GRAPH_DIR (there should be only one),
		then calculate the average PID between all pairs of seqs in <fasta_filename>.

		This is useful for calculating, say, the avg PID for Actinobacteria TPPs.
	"""
	pid_graph_filename = fnmatch.filter(os.listdir(PID_GRAPH_DIR),"*{0}*.sto.pid_graph.pickle".format(family))[0]
	pid_graph = read_gpickle(PID_GRAPH_DIR+'/'+pid_graph_filename)
	ids = []
	for line in os.popen("grep \">\" {0}".format(fasta_filename)):
		line = line.strip()[1:]
		db_id = line[(line.rfind('_')+1):]
		info = get_ncRNA_info(db_id)
		ids.append(gen_accID(info['acc'],info['start'],info['end'],info['strand'],info['acc_version']))

	g = pid_graph.subgraph(ids)
	pids = map(lambda e: e[2],g.edges_iter())
	return sum(pids)/len(pids)

def avg_pid_by_family_by_dir():
	rex = re.compile('\S+withDup_(\S+).fna')
	for filename in fnmatch.filter(os.listdir('.'),'*.fna'):
		try:
			family = rex.match(filename).group(1)
			pid = avg_pid_by_family(filename,family)
			print family,pid
		except:
			print >> sys.stderr, "failed on file ",filename

if __name__=="__main__":
	avg_pid_by_family_by_dir()

