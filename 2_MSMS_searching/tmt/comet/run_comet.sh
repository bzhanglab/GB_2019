#!/bin/sh

#SBATCH --job-name=comet
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --cpus-per-task=8
#SBATCH --time=100:00:00
#SBATCH --mem=60000M
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -e -x

batch_ID=${1}

### Initial compute node file system
if [ -d "/tmp/${batch_ID}" ];then
        rm -r "/tmp/${batch_ID}"
        mkdir "/tmp/${batch_ID}"
else
        mkdir "/tmp/${batch_ID}"
fi
if [ -d "/tmp/${batch_ID}/1st_step" ];then
        rm -r "/tmp/${batch_ID}/1st_step"
        mkdir "/tmp/${batch_ID}/1st_step"
else
        mkdir "/tmp/${batch_ID}/1st_step"
fi
if [ -d "/tmp/${batch_ID}/2nd_step" ];then
        rm -r "/tmp/${batch_ID}/2nd_step"
        mkdir "/tmp/${batch_ID}/2nd_step"
else
        mkdir "/tmp/${batch_ID}/2nd_step"
fi
if [ -d "/tmp/${batch_ID}/results" ];then
        rm -r "/tmp/${batch_ID}/results"
        mkdir "/tmp/${batch_ID}/results"
else
        mkdir "/tmp/${batch_ID}/results"
fi
if [ -d "/tmp/${batch_ID}/output" ];then
        rm -r "/tmp/${batch_ID}/output"
        mkdir "/tmp/${batch_ID}/output"
else
        mkdir "/tmp/${batch_ID}/output"
fi
if [ -d "/tmp/${batch_ID}/output/target_decoy" ];then
        rm -r "/tmp/${batch_ID}/output/target_decoy"
        mkdir "/tmp/${batch_ID}/output/target_decoy"
else
        mkdir "/tmp/${batch_ID}/output/target_decoy"
fi
if [ -d "/tmp/${batch_ID}/db" ];then
        rm -r "/tmp/${batch_ID}/db"
        mkdir "/tmp/${batch_ID}/db"
else
        mkdir "/tmp/${batch_ID}/db"
fi
if [ -d "/tmp/${batch_ID}/mgf" ];then
        rm -r "/tmp/${batch_ID}/mgf"
        mkdir "/tmp/${batch_ID}/mgf"
else
        mkdir "/tmp/${batch_ID}/mgf"
fi

### Download database file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/custom_database/${batch_ID}/ /tmp/${batch_ID}/db --recursive
### Download mgf file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/${batch_ID}.mgf /tmp/${batch_ID}/mgf
### Download mgf index file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/index_files/${batch_ID}.mgf_index.txt /tmp/${batch_ID}/mgf

### Comet search
/home/kail/software/comet/comet.2018014.linux.exe -P/home/kail/projects/2019_04_genome_biology/conf/scripts/comet/parameters/${batch_ID}_decoy.params -N/tmp/${batch_ID}/output/${batch_ID}_target_decoy /tmp/${batch_ID}/mgf/${batch_ID}.mgf
aws s3 cp /tmp/${batch_ID}/output s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/comet/${batch_ID}/ --recursive
mv /tmp/${batch_ID}/output/${batch_ID}_target_decoy.txt /tmp/${batch_ID}/results/
### Comet 1st step search
/home/kail/software/comet/comet.2018014.linux.exe -P/home/kail/projects/2019_04_genome_biology/conf/scripts/comet/parameters/${batch_ID}_1st.params -N/tmp/${batch_ID}/output/${batch_ID}_1st /tmp/${batch_ID}/mgf/${batch_ID}.mgf
aws s3 cp /tmp/${batch_ID}/output s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/comet/${batch_ID}/ --recursive
mv /tmp/${batch_ID}/output/${batch_ID}_1st.txt /tmp/${batch_ID}/results/

### Format comet results for pga
./format_comet_results.sh /tmp/${batch_ID}/results/${batch_ID}_target_decoy.txt
./format_comet_results.sh /tmp/${batch_ID}/results/${batch_ID}_1st.txt

### Caculate fdr
mkdir /tmp/${batch_ID}/results/global_fdr
mkdir /tmp/${batch_ID}/results/separate_fdr
Rscript pga_fdr_comet.R /tmp/${batch_ID}/results/${batch_ID}_target_decoy.txt_for_pga.txt /tmp/${batch_ID}/results/ "XXX_" /tmp/${batch_ID}/db/${batch_ID}_target_decoy.fasta
aws s3 cp /tmp/${batch_ID}/results/global_fdr s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/comet/${batch_ID}/global_fdr/ --recursive
aws s3 cp /tmp/${batch_ID}/results/separate_fdr s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/comet/${batch_ID}/separate_fdr/ --recursive

### Caculate 1st fdr and generate 2nd mgf
Rscript pga_fdr_comet_1st.R /tmp/${batch_ID}/results/${batch_ID}_1st.txt_for_pga.txt /tmp/${batch_ID}/1st_step/ "XXX_" /tmp/${batch_ID}/db/${batch_ID}_only_normal.fasta
aws s3 cp /tmp/${batch_ID}/1st_step/pga-peptideSummary.txt s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/comet/${batch_ID}/1st_step/
Rscript get_2st_mgf.R /tmp/${batch_ID}/1st_step/pga-peptideSummary.txt /tmp/${batch_ID}/mgf/${batch_ID}.mgf_index.txt /tmp/${batch_ID}/2nd_step/1st_step_spectrum_title.txt
java -Xmx60g -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar getSecondMGF /tmp/${batch_ID}/mgf/${batch_ID}.mgf /tmp/${batch_ID}/2nd_step/1st_step_spectrum_title.txt /tmp/${batch_ID}/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${batch_ID}/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/comet/${batch_ID}/2nd_step/ --recursive

### Run 2nd step identification
/home/kail/software/comet/comet.2018014.linux.exe -P/home/kail/projects/2019_04_genome_biology/conf/scripts/comet/parameters/${batch_ID}_2nd.params -N/tmp/${batch_ID}/2nd_step/${batch_ID}_2nd /tmp/${batch_ID}/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${batch_ID}/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/old_parameter_results/proteome_results/comet/${batch_ID}/2nd_step/ --recursive

### Clear node
rm -r /tmp/${batch_ID} 
