#!/bin/sh

set -e -x 

software=msgf+

table=/data/kail/2019_04_genome_biology/conf/data_table/prospective_colon/tmt_experiments.txt
sed 1d $table | while read -r tmtID
do	
	### Global data
	Rscript psm_fraction.R /data3/2019_04_genome_biology/${software}/${tmtID}/global_fdr/pga-peptideSummary.txt /data3/2019_04_genome_biology/mgf_files/index_files/${tmtID}.mgf_index.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_global_fdr.txt ${tmtID}
	### Sepatate data
	Rscript psm_fraction.R /data3/2019_04_genome_biology/${software}/${tmtID}/separate_fdr/pga-peptideSummary.txt /data3/2019_04_genome_biology/mgf_files/index_files/${tmtID}.mgf_index.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_separate_fdr.txt ${tmtID}
	### Two step data
	Rscript psm_fraction.R /data3/2019_04_genome_biology/two_step_search/${software}/2nd_step/${tmtID}/pga-peptideSummary.txt /data3/2019_04_genome_biology/two_step_search/${software}/2nd_step/mgf_index/${tmtID}_2nd.mgf_index.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_two_step.txt ${tmtID}
	Rscript combine_psm_pepquery.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_global_fdr.txt /data3/2019_04_genome_biology/pepquery_results/${tmtID}_pepquery.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_global_fdr_with_pepquery.txt 
	Rscript combine_psm_pepquery.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_separate_fdr.txt /data3/2019_04_genome_biology/pepquery_results/${tmtID}_pepquery.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_separate_fdr_with_pepquery.txt
	Rscript combine_psm_pepquery_2_step.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_two_step.txt /data3/2019_04_genome_biology/pepquery_results/${tmtID}_pepquery.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/${software}/${tmtID}_two_step_with_pepquery.txt
done
	

