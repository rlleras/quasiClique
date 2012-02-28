#!/bin/bash

# #######################################################
# Similar to run_pscore_pfold.sh except that
# I'm using the .tree from pfold as the starting
# tree for running RAxML to get a ML tree
# then finally running pscore on that
# 
# NOTE: requires pscore_pfold to have been run ALREADY
# (uses the .fasta and .tree file!)
# 
# #######################################################

baseDir=/external3/home/etseng/software_downloads/
cmf="$baseDir/CMfinder_03"

FAStoPHY="/external3/home/etseng/lib/python/convert/fasta_to_phyaln.py"
RAxML="$baseDir/RAxML-7.0.4/raxmlHPC"
PSCORE="$cmf/c/pscorez/posterior"

dirname=$1

cd $dirname
for FILE in `ls *.motif*.h[1-9]_[1-9]`
do
	basename=$FILE
	# step 1: converting .fasta to .phy_aln
	python $FAStoPHY -f $basename.fasta -o $basename.phy_aln
	# step 2: use .tree from pfold as initial tree for RAxML
	#         the result tree is RAxML_result.$basename
	$RAxML -m GTRGAMMA -s $basename.phy_aln -n $basename -t $basename.tree
	# step 3: run pscore on the RAxML tree
	$PSCORE --partition -t RAxML_result.$basename -s $cmf/data/single_model -p $cmf/data/pair_model -g $cmf/data/scfg $basename > $basename.pscoreout_RAxML
done
