import os
import sys
import csv
import glob
from cPickle import *
from collections import defaultdict

MOTIF_DIR = 'motifs/'

STORE_DIR = 'MotifScore_vs_CMhit/'

REAL_NCRNA_SS_CONS_FILE = '/home/etseng/silo/UW_Larry/data/2009-08-19-sto/ALLncRNA.SS_cons.txt'

CM_DB_DIR = '/home/etseng/silo/UW_Larry/data/CMscanDB_100X/' # '/external3/home/etseng/CMscanDB_100X'
CM_DB_SUMMARY = os.path.join(CM_DB_DIR, 'SUMMARY')
# read the summary file to get the family sizes
db_summary = {}
for obj in csv.DictReader(open(CM_DB_SUMMARY), delimiter='\t'):
		db_summary[obj['FAMILY']] = {'TRUE': int(obj['TRUE']), 'CONTROLS': int(obj['CONTROLS'])}

# read the ss_cons file to get (Zasha's annotated) SS_cons for each family
for line in open(REAL_NCRNA_SS_CONS_FILE):
	family, ss_cons = line.strip().split('\t')
	if family not in db_summary: continue
	db_summary[family]['SS_CONS'] = ss_cons

with open('description.pickle') as f: dir_info = load(f)


