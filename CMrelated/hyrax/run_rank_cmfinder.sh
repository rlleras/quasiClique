#!/bin/bash

baseDir=/external3/home/etseng/software_downloads/
script="$baseDir/NewMotifScore.pl"
cmf="$baseDir/CMfinder_03"

HYRAXCLUSTER='uplow'

for FILE in `ls -d Firm_HC_*`
do
	DIR=`basename $FILE`
	echo "cd $DIR"
	echo "srun -p $HYRAXCLUSTER -N 1 rank_cmfinder.pl -w -rank \"${DIR}.fna.[1-9].motif*\" ${DIR}.fna.summary &"
	echo "cd .."
done
