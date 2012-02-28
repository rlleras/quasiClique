import os,re,sys,fnmatch
from collections import defaultdict
from cPickle import *

def merge_pickles(pickle_list, output_filename, delete_merged):
	QQQ = []
	for file in pickle_list:
		with open(file) as handle:
			QQQ += load(handle)
		if delete_merged:
			os.system("rm " + file)

	with open(output_filename, 'w') as handle:
		dump(QQQ, handle)

if __name__ == "__main__":
	from optparse import OptionParser
	parser = OptionParser()
	parser.add_option("-p", "--prefix", help="pickle prefix so pickle is form <prefix>_<machine>_<index>.pickle")
	parser.add_option("-d", "--deleteMerged", help="delete the merged pickles", action="store_true", default=False)

	(options, args) = parser.parse_args()

	rex = re.compile(options.prefix+'_(\S+)_(\d+).pickle')
	machines = defaultdict(lambda: [])

	for file in fnmatch.filter( os.listdir('.'), options.prefix+'*.pickle' ):
		m = rex.match(file)
		machine,index = m.group(1),m.group(2)
		# ignore empty pickles
		if os.stat(file).st_size > 0:
			machines[machine].append( file )

	print("------ starting to merge... ------")
	print("Deletes merged pickles? {0}".format( options.deleteMerged ))

	for machine, pickle_list in machines.iteritems():
		pickle_list.sort()
		output_filename = "{prefix}_{machine}_ALL.pickle".format(prefix=options.prefix,\
				machine=machine)
		print("Merging {0}\n to {1}".format("\n".join(pickle_list), output_filename))
		raw_input("Go head? press any to continue")
		merge_pickles(pickle_list, output_filename, options.deleteMerged)

