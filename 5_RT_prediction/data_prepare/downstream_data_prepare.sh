#!/bin/sh
set -e -x

software=msgf+

table=/data/kail/2019_04_genome_biology/conf/data_table/prospective_colon/tmt_experiments.txt
sed 1d $table | while read -r tmtID
do
	mkdir /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/${software}/${tmtID}
	Rscript psm_fraction.R /data3/2019_04_genome_biology/${software}/${tmtID}/global_fdr/pga-peptideSummary.txt /data/kail/2019_04_genome_biology/data/prospective_colon/mgf_files/index_files/${tmtID}.mgf_index.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/${software}/${tmtID}/${tmtID}_fraction
	python prepare_train_data.py /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/${software}/${tmtID}/${tmtID}_fraction_normal_psm.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/${software}/${tmtID}/

done


