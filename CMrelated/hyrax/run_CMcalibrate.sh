#!/bin/bash

# #####################################################################
# new Infernal (> 1.0) doesn't accept CMFinder's old cm format
# so must first run cmbuild to get .cm files
# #####################################################################

CMD='cmcalibrate'
HYRAXCLUSTER='pubnorm'

for FILE in `ls */*.motif*.h[1-9]_[1-9]*.cm`
do
	dirname=`dirname $FILE`
	basename=`basename $FILE`
	echo "srun -N 1 -p $HYRAXCLUSTER $CMD $FILE&"
done
