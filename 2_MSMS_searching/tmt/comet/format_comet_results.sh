#!/bin/sh

set -e -x
comet_result=${1}

first=$(head -n 1 $comet_result)

if [[ $first == *"CometVersion"* ]]
then
	sed -i '1d' $comet_result
	sed "s/modifications/modifications\tna/g" $comet_result > ${comet_result}.tmp
	mv ${comet_result}.tmp $comet_result
fi

Rscript format_comet_results.R $comet_result ${comet_result}_for_pga.txt
