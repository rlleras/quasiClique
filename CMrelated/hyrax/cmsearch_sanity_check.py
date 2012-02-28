import os
import sys
import glob

SHOULD_SEE = [\
'command:',\
'fil-T-hmm',\
'target name',\
'Post-search info for',\
'run time']

def cmsearched_sanity_check(filename):
	with open(filename) as f:
		chunk = f.read()

	test = [chunk.find(s) for s in SHOULD_SEE]
#	for s in SHOULD_SEE: print s, chunk.find(s)
#	raw_input()
	if all((test[i] < test[i+1] and test[i]>=0) for i in xrange(len(SHOULD_SEE)-1)):
		pass
	else:
		print("{0} DID NOT PASS SANITY CHECK!!!!!!!!!!".format(filename))

if __name__ == "__main__":
	for d in os.listdir('.'):
		for file in glob.iglob("{d}/*.motif.*h[1-9]_[1-9]".format(d=d)):
			cmsearched_sanity_check(file+'.cmsearched')

