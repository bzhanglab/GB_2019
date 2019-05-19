#!/bin/sh

set -e -x 

table=/data/kail/2019_04_genome_biology/conf/data_table/prospective_colon/tmt_experiments.txt

batch_ID=""
sed 1d $table | while read -r tmtID
do
        batch_ID=${tmtID}
        echo ${batch_ID}
	mkdir /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/rt_cobination/${batch_ID}
	for file in /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/msgf+/${batch_ID}/*_normal.txt
	do
		comet=`echo "${file/msgf+/comet}"`
		xtandem=`echo "${file/msgf+/xtandem}"`
		output=`echo "${file/msgf+/rt_cobination}"`
		Rscript rt_normal_combination.R $file $xtandem $comet $output
	done
done
