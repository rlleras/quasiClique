import os,re,sys
from miscUtils import line_count
from timeit import Timer
from getopt import getopt

if len(sys.argv) < 4:
	print >> sys.stderr, "\nUSAGE: python {0} <query_fasta> <dbname> <output_dir> <extra_output_prefix>\n".format(sys.argv[0])
	sys.exit(-1)

query_fasta_filename = sys.argv[1]
dbname = sys.argv[2]
output_dir = sys.argv[3]
extra_prefix = ''
if len(sys.argv) > 4: extra_prefix=sys.argv[4]
if len(sys.argv) > 5: dynamic_args=sys.argv[5:]

M = "5"
N = "4"
Q = "10"
R = "6"
W = "3"
E = "10"
optlist, args = getopt(dynamic_args,'M:N:Q:R:W:E:')
for o,v in optlist:
	if o=='-M': M=v
	elif o=='-N': N=v
	elif o=='-Q': Q=v
	elif o=='-R': R=v
	elif o=='-W': W=v
	elif o=='E': E=v

output_prefix = "M{0}N{1}Q{2}R{3}W{4}E{5}".format(M,N,Q,R,W,E)
output_cmd = "-M {0} -N {1} -Q {2} -R {3} -W {4} -E {5}".format(M,N,Q,R,W,E)
if len(sys.argv) > 4: output_prefix = output_prefix+'.'+extra_prefix

if not os.path.exists(dbname + '.xnd'):
	os.system("xdformat -n -o {0} {1}".format(dbname,dbname))

total_runtime = 0

with open("{2}/{0}.{1}.WUblast.timing".format(os.path.basename(query_fasta_filename),output_prefix,output_dir),'w') as timing_f:
	for i in xrange(line_count(query_fasta_filename,'fasta')-1):
	#       os.system("blastn {0} {1} -M 5 -N -4 -Q 10 -R 6 -W 7 -E 10 -noseqs -cpus 4 -mformat 2 -o {2}.M5N4Q10R6W7E10.WUblast_part{3} -qrecmin {4} -qrecmax {5} -dbrecmin {6}".format(dbname,query_fasta_filename,query_fasta_filename,i,i+1,i+1,i+2))
		cmd = "blastn {0} {1} {2} -noseqs -cpus 4 -mformat 2 -o {8}/{3}.{4}.WUblast_part{5} -qrecmin {6} -qrecmax {7}".format(dbname,query_fasta_filename,output_cmd,os.path.basename(query_fasta_filename),output_prefix,i,i+1,i+1,output_dir)
#		if proportinateDB is not None:
#			size = int(os.popen("grep -c \">\" " + query_fasta_filename).readline().strip())
#			cmd += " -dbrecmax " + str(size*proportinateDB)
		t = Timer("os.system('{0}')".format(cmd),'import os')
		runtime = t.timeit(1)
		total_runtime += runtime
		timing_f.write("part{0}\t{1}\n".format(i,runtime))
		sys.stdout.flush()
	timing_f.write("TOTAL\t"+str(total_runtime)+"\n")
	tmp = os.path.basename(query_fasta_filename)+'.'+output_prefix
	os.system("python concat_blast_files.py {0} {1} {2}".format(output_dir,tmp+'.WUblast_part',output_dir+'/'+tmp+'.WUblast',True))
