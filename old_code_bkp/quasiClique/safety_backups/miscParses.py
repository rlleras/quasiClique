import os,re,sys

#id_rex = re.compile('(\S+)\\/([\-]*\d+)-([\-]*\d+)')

def adjust_coord((strand,start,end),loc_start,loc_end):
	if start >= end:
		raise Exception, "start must be < end! (adjust_coord)"
	if strand == +1:
		if loc_start < loc_end:
			return (+1,start+loc_start-1,start+loc_end-1)
		else:
			return (-1,start+loc_end-1,start+loc_start-1)
	else:
		if loc_start < loc_end:
			return (-1,end-loc_end+1,end-loc_start+1)
		else:
			return (+1,end-loc_start+1,end-loc_end+1)

def findOriginalNode(acc,strand,start,end,table_prefix,cursor):
        mysql_query = 'select id,start,end from ' + table_prefix + '_node_parsed where acc="' + acc + '" and strand=' + str(strand) + ' and does_overlap(' + str(start) + ',' + str(end) + ',start,end) order by calc_overlap(' + str(start) + ',' + str(end) + ',start,end) desc limit 1'
        cursor.execute(mysql_query)
        data = cursor.fetchone()
        return [int(data[0]),int(data[1]),int(data[2])]

def uniquify1(seq):
	set = {} 
        map(set.__setitem__, seq, []) 
        return set.keys()

def calc_overlap(a1,a2,b1,b2):
        if a1 <= b1:
                return min(b2,a2)-b1+1
        else: # b1 < a1
                return min(b2,a2)-a1+1

def calc_interval_overlap(interval1,interval2):
	return calc_overlap(interval1.lower_bound,interval1.upper_bound,interval2.lower_bound,interval2.upper_bound)

def gen_accID(acc,start,end,strand,version=None):
	if version is not None and version!='':
		acc = acc+'.'+str(version)
	if int(strand)==+1:
		return "{0}/{1}-{2}".format(acc,start,end)
	else:
		return "{0}/{1}-{2}".format(acc,end,start)

def split_acc_version(acc):
	i = acc.find('.')
	if i > 0:
		return (acc[:i],acc[(i+1):])
	return (acc,None)
	
def remove_acc_version(acc):
	i = acc.find('.')
	if i > 0:
		return acc[:i]
	return acc	

def parsed_accID_nostrand(id):
	"""
		Unlike parsed_accID, this does not take the loc_start,loc_end argument
		It also doesn't return the strand. Instead it simple returns
		(acc,start,end) s.t. id is acc/start-end, without caring whether start < end
	"""
	i = id.find('/')
	j = id.find('-')
	return (id[:i],int(id[(i+1):j]),int(id[(j+1):]))

def parsed_accID(id,version_split=False,loc_start=None,loc_end=None):
	"""
	pre: loc_start and loc_end should either both be None or not None.
	id should be in form accession/global_start-global_end,
	returns (acc,strand,start,end)
	"""
	# re too slow? let's try this instead
	i = id.find('/')
	j = id.find('-')
	acc = id[:i]
	if version_split:
		acc = split_acc_version(acc)
	a = int(id[(i+1):j])
	b = int(id[(j+1):])		
	if loc_start is None and loc_end is None:
		if a < b:
			return (acc,+1,a,b)
		else: 
			return (acc,-1,b,a)
	else:
	        strand,start,end = adjust_coord((1,a,b),loc_start,loc_end)
		return (acc,strand,start,end)

def parseWUBLASTline(line):
	batch = line.rstrip().split('\t')
	result = {'id1':batch[0],
	          'id2':batch[1],
                  'e':float(batch[2]),
                  'sprime':float(batch[4]),
		  's':float(batch[5]),
                  'pcident':float(batch[10]),
                  'pcpos':float(batch[11]),
                  'strand1':int(batch[16]),
                  'start1':int(batch[17]),
                  'end1':int(batch[18]),
                  'strand2':int(batch[19]),
                  'start2':int(batch[20]),
                  'end2':int(batch[21])}
	return result

def parseNCBIBLASTline(line):
	batch = line.rstrip().split('\t')
	return {'id1':batch[0],
		'id2':batch[1],
		'pcident':float(batch[2]),
		'length':int(batch[3]),
		'pcmismatch':float(batch[4]),
		'gaps':int(batch[5]),
		'start1':int(batch[6]),
		'end1':int(batch[7]),
		'start2':int(batch[8]),
		'end2':int(batch[9]),
		'e':float(batch[10]),
		'sprime':float(batch[11])}

# pre: the main cluster id should be the same
def matchClusterPrefix(sub_prefix1,sub_prefix2):
	if sub_prefix1==sub_prefix2:
		return len(sub_prefix1)
	for i in range(0,min(len(sub_prefix1),len(sub_prefix2))):
		if sub_prefix1[i]!=sub_prefix2[i]:
			return i

# input: list of cluster ids
# output: calculated averaged pairwise distance
def calculateFamilyDistance(C):
	if len(C)==0 or len(C)==1:
		return float('Inf')

	new_C = []
	# first, parse into (main_id,sub_prefix)
	for c in C:
		new_C.append((c[:c.find('.')],c[(c.find(':')+1):]))
	
	# make up a maximum distance for members with different main_id
	# NOTE!!! MAX_DIST set to 70!!!
	MAX_DIST = 70.

	(dist_acc,num_acc) = (0,0)
	for i in range(len(new_C)):
		for j in range(i+1,len(new_C)):
			num_acc += 1
			if new_C[i][0]!=new_C[j][0]: # different main_id
				dist_acc += MAX_DIST
			else:
				match_prefix_length = matchClusterPrefix(new_C[i][1],new_C[j][1])
				dist_acc += min(MAX_DIST, len(new_C[i][1]) + len(new_C[j][1]) - 2 * match_prefix_length)
		
	return dist_acc * 1. / num_acc

# input: C is a list of names
def findClosestSubcluster(query,C):
	# first, only pick out C members that have the same main_id as query
	main_id = query[:query.find('.')]
	newC = []
	for c in C:
		if c[:c.rfind('.')]==main_id:
			newC.append(c)
	C = newC # swap	
	if len(C)==1:
		return (None,None)

	C.sort()
	# then, find query position in C
	pos = C.index(query)
	# the closest one is either pos-1 or pos+1
	# tend to special boundary case first
	if pos==0:
		return (matchClusterPrefix(query,C[1]),C[1])
	elif pos==len(C)-1:
		return (matchClusterPrefix(query,C[len(C)-2]),C[len(C)-2])
	else:
		prevD = matchClusterPrefix(query,C[pos-1])
		nextD = matchClusterPrefix(query,C[pos+1])
		if prevD==nextD:
			#print 'Tie between ' + C[pos-1] + ' and ' + C[pos+1]
			return (None,None)
		elif prevD > nextD:
			return (prevD,C[pos-1])
		else:
			return (nextD,C[pos+1])


def msgRiboHit(assoc_type):
	return_msg = ''
	for type,unique_hits in assoc_type.iteritems():
		return_msg += ',' + str(len(unique_hits)) + ' ' + type

	return return_msg[1:]


def read_motif(motif_filename):
	"""
		Reads a motif file (probably outputted from CMFinder)
		which is in Stockholm format (but doesn't check for the
		silly STOCKHOLM header)

		Takes the SS_cons consensus motif structure, and reformats
		it so it only consists of ( ) .
		
		Returns (consensus_sequence,consensus_structure) with all
		gaps removed
	"""
	from Bio import SeqIO
	valid_left  = ('(','<','{','[')
	valid_right = (')','>','}',']')
	
	final_rf = ''
	final_ss = ''
	
	try:
		rec = SeqIO.parse(open(motif_filename),'stockholm')
	except:
		return (None,None)
		
	for i,x in enumerate(rec.rf):
		if x <> '.':
			final_rf += x.upper()
			if rec.ss_cons[i] in valid_left:
				final_ss += '('
			elif rec.ss_cons[i] in valid_right:
				final_ss += ')'
			else:
				final_ss += '.'
	return (final_rf,final_ss)
			
