import os,re,sys
from Bio import SeqIO
from random import random
from bisect import *

def no_ambiguous_code(s):
	for i in xrange(len(s)):
		if s[i] not in ['A','U','T','C','G']:
			return False
	return True

ribo_filename = sys.argv[1] # should be in fasta format
igr_dirname = sys.argv[2]
dice_freq = float(sys.argv[3])
min_len = 100
max_len = 5000

ribo_acc_set = {}
rex = re.compile('(\S+)\\/(\d+)-(\d+)')
for r in SeqIO.parse(open(ribo_filename),'fasta'):
	try:
		m = rex.match(r.id)
		acc = m.group(1)
		i = acc.find('.')
		if i > 0: acc = acc[:i]
		if acc not in ribo_acc_set: ribo_acc_set[acc] = []
		insort(ribo_acc_set[acc],(int(m.group(2)),int(m.group(3))))
	except:
		print >> sys.stderr, "failed extracting acc info for %s" % r.id

for file in os.listdir(igr_dirname):
	#acc = file[:file.find('.')]
	for r in SeqIO.parse(open(igr_dirname+'/'+file),'fasta'):
		m = rex.match(r.id)
		acc = m.group(1)
		i = acc.find('.')
		if i > 0: acc = acc[:i]
		
		ok_random = False
		if acc in ribo_acc_set:
			if no_ambiguous_code(r.seq.tostring()):
				start = int(m.group(2))
				end = int(m.group(3))
				if min_len <= end-start <= max_len:
					i = bisect(ribo_acc_set[acc],(start,end))
					if i == 0:
						ok_random = end < ribo_acc_set[acc][0][0]
					elif i == len(ribo_acc_set[acc])-1:
						ok_random = ribo_acc_set[acc][-1][1] < start
					elif i <= len(ribo_acc_set[acc])-2:
						ok_random = ribo_acc_set[acc][i-1][1] < start and end < ribo_acc_set[acc][i][0]
					else:
						ok_random = ribo_acc_set[acc][i-1][1] < start
		else:
			ok_random = True
		if ok_random:
			print ">%s\n%s" % (r.id,r.seq.tostring())

