#!/bin/sh

set -e -x

software=msgf+

table=/data/kail/2019_04_genome_biology/conf/data_table/prospective_colon/tmt_experiments.txt
sed 1d $table | while read -r tmtID
do
	Rscript prepare_pepquery_input.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/${software}/${tmtID}/${tmtID}_fraction_variant_psm.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/${software}/pepquery_input/${tmtID}_global_variant_peptide.txt
done
