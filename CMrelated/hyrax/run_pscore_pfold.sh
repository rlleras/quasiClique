#!/bin/bash

baseDir=/external3/home/etseng/software_downloads/
script="$baseDir/NewMotifScore.pl"
cmf="$baseDir/CMfinder_03"
pfold="$baseDir/pfold/"
pscoreRate="$baseDir/pfold/pscore.rate"

HYRAXCLUSTER='uplow'

for FILE in `ls */*.motif*.h[1-9]_[1-9]`
do
	dirname=`dirname $FILE`
	basename=`basename $FILE`
	echo "srun -N 1 -p HYRAXCLUSTER perl $script -tempBase $dirname/$basename -cmfinderDir $cmf -pfoldPath $pfold $dirname/$basename pscore-pfold-sto $pscoreRate&"
done
