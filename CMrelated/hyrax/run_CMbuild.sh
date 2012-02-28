#!/bin/bash

# ##################################################################
# this only runs cmbuild (which is fast)
# if want to get E-values use cmcalibrate <cmfile> (VERY SLOW)
# ##################################################################

CMD='cmbuild'
HYRAXCLUSTER='uplow'

for FILE in `ls */*.motif*.h[1-9]_[1-9]`
do
	dirname=`dirname $FILE`
	basename=`basename $FILE`
	echo "srun -N 1 -p $HYRAXCLUSTER $CMD $FILE.cm $FILE&"
done
