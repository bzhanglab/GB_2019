#!/bin/sh

set -e -x

software=msgf+
level=peptide_level

Rscript combine_psm_fraction.R /data3/2019_04_genome_biology/${software}/combination_fdr/peptide_level/global_fdr/pga-peptideSummary.txt /data3/2019_04_genome_biology/mgf_files/index_files/all_index_information.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/combine_fdr_results/combine_global_fdr_${software}.txt

Rscript combine_psm_fraction.R /data3/2019_04_genome_biology/two_step_search/${software}/2nd_step/combination_fdr/${level}/global_fdr/pga-peptideSummary.txt /data3/2019_04_genome_biology/two_step_search/${software}/2nd_step/mgf_index/all_index_information.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/${level}/combine_fdr_results/two_step_${software}.txt

Rscript combine_psm_fraction.R /data3/2019_04_genome_biology/${software}/combination_fdr/peptide_level/separate_fdr/pga-peptideSummary.txt /data3/2019_04_genome_biology/mgf_files/index_files/all_index_information.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/combine_fdr_results/combine_separate_fdr_${software}.txt

Rscript combine_combine_psm_pepquery.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/combine_fdr_results/combine_fdr_peptide_level_${software}_global_fdr.txt /data3/2019_04_genome_biology/pepquery_results/all_tmt_pepquery.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/combine_fdr_results/combine_global_fdr_with_pepquery_${software}.txt
Rscript combine_combine_psm_pepquery.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/combine_fdr_results/combine_fdr_peptide_level_${software}_separate_fdr.txt /data3/2019_04_genome_biology/pepquery_results/all_tmt_pepquery.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/peptide_level/combine_fdr_results/combine_separate_fdr_with_pepquery_${software}.txt
Rscript combine_combine_psm_pepquery.R /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/${level}/combine_fdr_results/two_step_${software}.txt /data3/2019_04_genome_biology/pepquery_results/all_tmt_pepquery.txt /data/kail/2019_04_genome_biology/data/prospective_colon/downstream_analysis/new_format/${level}/combine_fdr_results/combine_fdr_${level}_${software}_two_step.txt


