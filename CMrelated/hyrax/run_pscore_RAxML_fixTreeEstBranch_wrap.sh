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

for DIRNAME in `ls -d */`
do
	echo "srun -p publow -N 1 bash run_pscore_RAxML_fixTreeEstBranch_single.sh $DIRNAME&"
done
