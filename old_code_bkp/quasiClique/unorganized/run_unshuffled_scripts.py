import os,re,sys,fnmatch
from rename_shuffled_seqs import squid_shuffle_seq

def test1(input_dirname,output_dirname,X=60,dynamic_param=''):
	"""
		Run each family independently, where
			query - ncRNAs from this family
			DB - query + each query dishuffled <X> times

		The resulting files will have a pattern of _plus_self<X>_modishuffled
	"""
	print >> sys.stderr, X
	for filename in fnmatch.filter(os.listdir(input_dirname),'*.fna'):
		query_filename = input_dirname+'/'+filename
		print >> sys.stderr, query_filename,'.......'
		dishuffled_filename = output_dirname+'/'+filename[:-4]+'_self'+str(X)+'_modishuffled.fna'
		DB_name = filename[:-4]+'_plus_self'+str(X)+'_modishuffled'
		#print("shuffle -d -n {2} -o {0} {1}".format(dishuffled_filename,query_filename,X))
		squid_shuffle_seq(query_filename,dishuffled_filename,X)
		os.system("xdformat -n -o {0} {1}".format(DB_name,query_filename))
		os.system("xdformat -n -a {0} {1}".format(DB_name,dishuffled_filename))
		os.system("python run_unshuffled_blast_dynamic_param.py {0} {1} {2} plus_self{3}_modishuffled {4}".format(query_filename,DB_name,output_dirname,X,dynamic_param))

def test2(input_dirname,dishuffled_filename,output_dirname,X,dynamic_param):
	"""
		Run each family independently, where
			query - ncRNAs from this family
			DB - query + propotional chunk of dishuffled seq.
				    (make sure this dishuffled seq is ONE line ID followed by ONE line seq!)

		The resulting files will have a pattern of _plus_prop<X>_modishuffled
	"""
	for filename in fnmatch.filter(os.listdir(input_dirname),'*.fna'):
                query_filename = input_dirname+'/'+filename
                print >> sys.stderr, query_filename,'.......'
		query_size = int(os.popen("grep -c \">\" " + query_filename).readline().strip())
		DB_name = filename[:-4]+'_plus_prop'+str(X)+'_modishuffled'
		os.system("xdformat -n -o {0} {1}".format(DB_name,query_filename))
		tmp = os.tempnam()
		# the following line works ONLY if the fasta is one line for ID, one line for seq!!!
		os.system("head -n {0} {1} > {2}".format(query_size*X*2,dishuffled_filename,tmp))
		os.system("xdformat -n -a {0} {1}".format(DB_name,tmp))
		os.system("python run_unshuffled_blast_dynamic_param.py {0} {1} {2} plus_prop{3}_modishuffled {4}".format(query_filename,DB_name,output_dirname,X,dynamic_param))

def test3(input_dirname,random_filename,output_dirname):
	"""
		Run each family independently, where
			query - ncRNAs from this family
			DB - query + random_file
	"""
	for filename in fnmatch.filter(os.listdir(input_dirname),'*withDup_y*.fna'):
                query_filename = input_dirname+'/'+filename
                print >> sys.stderr, query_filename,'.......'
		DB_name = filename[:-4]+'_plus_random'
		os.system("xdformat -n -o {0} {1}".format(DB_name,query_filename))
		os.system("xdformat -n -a {0} {1}".format(DB_name,random_filename))
		os.system("python run_unshuffled_blast.py {0} {1} {2} plus_random".format(query_filename,DB_name,output_dirname))

def test4(input_dirname,random_filename,output_dirname,X=60,dynamic_param=''):
        for filename in fnmatch.filter(os.listdir(input_dirname),'*.fna'):
                query_filename = input_dirname+'/'+filename
                print >> sys.stderr, query_filename,'.......'
                query_size = int(os.popen("grep -c \">\" " + query_filename).readline().strip())
                DB_name = filename[:-4]+'_plus_prop'+str(X)+'_random_'
                if not os.path.exists(DB_name+".xnd"):
                        os.system("xdformat -n -o {0} {1}".format(DB_name,query_filename))
                        tmp = os.tempnam()
                        # the following line works ONLY if the fasta is one line for ID, one line for seq!!!
                        os.system("head -n {0} {1} > {2}".format(query_size*X*2,random_filename,tmp))
                        os.system("xdformat -n -a {0} {1}".format(DB_name,tmp))
                os.system("python run_unshuffled_blast_dynamic_param.py {0} {1} {2} plus_prop{3}_random {4}".format(query_filename,DB_name,output_dirname,X,dynamic_param))

def test5(input_dirname,dishuffled_filename,output_dirname,X=60,dynamic_param=''):
        for filename in fnmatch.filter(os.listdir(input_dirname),'*.fna'):
                query_filename = input_dirname+'/'+filename
                print >> sys.stderr, query_filename,'.......'
                query_size = int(os.popen("grep -c \">\" " + query_filename).readline().strip())
                DB_name = filename[:-4]+'_plus_prop'+str(X)+'_modishuffled_'
		if not os.path.exists(DB_name+".xnd"):
	                os.system("xdformat -n -o {0} {1}".format(DB_name,query_filename))
        	        tmp = os.tempnam()
	                # the following line works ONLY if the fasta is one line for ID, one line for seq!!!
	                os.system("head -n {0} {1} > {2}".format(query_size*X*2,dishuffled_filename,tmp))
	                os.system("xdformat -n -a {0} {1}".format(DB_name,tmp))
		real_q_filename = query_filename.replace('_padseqlen','')
                os.system("python run_unshuffled_blast_dynamic_param.py {0} {1} {2} plus_padseqlen_prop{3}_modishuffled {4}".format(real_q_filename,DB_name,output_dirname,X,dynamic_param))

if __name__=="__main__":
#	test1(sys.argv[1],sys.argv[2],int(sys.argv[3]),' '.join(sys.argv[4:]))
#	test2(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),' '.join(sys.argv[5:]))	
	test5(sys.argv[1],sys.argv[2],sys.argv[3],int(sys.argv[4]),' '.join(sys.argv[5:]))
#	test3(sys.argv[1],sys.argv[2],sys.argv[3])
