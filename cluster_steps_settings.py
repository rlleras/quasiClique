import os,re,sys,timeit,time,fnmatch,itertools
from cPickle import *
from bisect import *
from collections import defaultdict
# non-standard or self-written packages
from networkx import *
from networkx.readwrite import *
from miscParses import calc_overlap, parseWUBLASTline, parseNCBIBLASTline, adjust_coord, parsed_accID
from miscBio import FastaReader
from interval import *
from operator import itemgetter
from random import uniform

try:
	from miscMySQL import *
	CONN_FUNC = get_conn_Harwood
except ImportError:
	print >> sys.stderr, "not importing miscMySQL"

try:
	from psyco.classes import *		
except ImportError:
	print >> sys.stderr, "not importing psyco.classes"
	
CLIQUE_MIN_SIZE = 5
PERFECT_CLIQUE_MIN_SIZE = 3
PERFECT_MAXITR = 20
QUASI_GAMMA = 0.8
QUASI_MAXITR = 20
