import os,re,sys,miscMySQL,miscTaxonomy
from collections import defaultdict
from decimal import Context, localcontext

conn = miscMySQL.get_conn_ncRNA()
d_cur = miscMySQL.get_dict_cursor(conn)
connNCBI = miscMySQL.get_conn_NCBI()
curNCBI = connNCBI.cursor()

def tally():
	result = defaultdict(lambda: defaultdict(lambda: ''))
	d_cur.execute("select family,count(*) as total from Zasha_20081002_plus20071102_curated group by family")
	for r in d_cur.fetchall():
		result[r['family']]['total'] = r['total']
	d_cur.execute("select z.family as family,c.phylum as phylum, count(*) as count from Zasha_20081002_plus20071102_curated as z left join NCBI.cache_acc_to_tax as c on (z.acc=c.acc) group by z.family,c.phylum")
	for r in d_cur.fetchall():
		result[r['family']][r['phylum']] = r['count']

	d_cur.execute("select distinct phylum from NCBI.cache_acc_to_tax order by phylum")
	phyla = map(lambda x: x['phylum'],d_cur.fetchall())
	print("\tTOTAL\t"+"\t".join(phyla))
	for family,d in result.iteritems():
		print("{0}\t{1}".format(family,d['total'])),
		for p in phyla:
			print("\t"+str(d[p])),
		print ''

tally()
sys.exit(-1)

# here's what we're interested in:
# 1. Total sum of ncRNA nts per phylum
# 2. Average Ratio of ncRNA/IGR (rough estimate) per phylum
bin_total_ncRNA = defaultdict(lambda: 0)
bin_ratio = defaultdict(lambda: [])

d_cur.execute("select acc,sum(end-start+1) as s from Zasha_20081002_plus20071102_curated group by acc order by s desc")
for result in d_cur.fetchall():
	# will get None is this is not a complete genome in RefSeq
	tax = miscTaxonomy.get_phylum(result['acc'])
	if tax is None: continue
	ncRNA_length = int(result['s']) 
	if 'phylum' in tax:
		# get the IGR size
		try:
			curNCBI.execute("select length from RefSeq25_IGRs_m30s0_sizes where acc_w_version like '{0}%'".format(result['acc']))
			IGR_length = int(curNCBI.fetchone()[0])
		except:
			continue
		bin_total_ncRNA[tax['phylum']] += ncRNA_length
		bin_ratio[tax['phylum']].append(ncRNA_length*1./IGR_length)

print '============ total ncRNA nts (doncare if not in IGR) per phylum ============'
for k,v in bin_total_ncRNA.iteritems():
	print k,v
print '============ median ratio of ncRNA/IGR per phylum ============'
for k,v in bin_ratio.iteritems():
	v.sort()
	l = len(v)
	print k,"[{0: .4f},{1: .4f},{2: .4f}]".format(v[l/4],v[l/2],v[3*l/4])

d_cur.close()
conn.close()
curNCBI.close()
connNCBI.close()
