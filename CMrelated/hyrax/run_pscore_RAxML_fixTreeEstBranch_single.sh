#!/bin/bash

# #######################################################
# Similar to run_pscore_RAxML.sh but
# (1) uses the tree from pfold and does NOT alter it!
# (2) simply estimate branch length (-f e option)
#
# NOTE: requires pscore_pfold to have been run ALREADY
# (uses the .fasta and .tree file!)
#
# NOTE: to differentiate the output from running run_pscore_RAxML.sh
#       set output -n <dirname>/<basename>.fixTreeEstBranch
#
# #######################################################

baseDir=/external3/home/etseng/software_downloads/
cmf="$baseDir/CMfinder_03"

FAStoPHY="/external3/home/etseng/lib/python/convert/fasta_to_phyaln.py"
RAxML="$baseDir/RAxML-7.0.4/raxmlHPC"
PSCORE="$cmf/c/pscorez/posterior"

dirname=$1
d_b=`basename $dirname`

cd $dirname
for FILE in `ls ${d_b}.*.motif*.h[1-9]_[1-9]`
do
	basename=$FILE
	echo $basename.........
	# step 1: converting .fasta to .phy_aln
	python $FAStoPHY -f $basename.fasta -o $basename.phy_aln
	# step 2: use .tree from pfold as initial tree for RAxML
	#         the result tree is RAxML_result.$basename
	$RAxML -f e -m GTRGAMMA -s $basename.phy_aln -n $basename.fixTreeEstBranch -t $basename.tree
	# step 3: run pscore on the RAxML tree
	$PSCORE --partition -t RAxML_result.$basename.fixTreeEstBranch -s $cmf/data/single_model -p $cmf/data/pair_model -g $cmf/data/scfg $basename > $basename.fixTreeEstBranch.pscoreout_RAxML
done
