#!/bin/sh

#SBATCH --job-name=pepquery
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=55000M
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -e -x

sample_ID=${1}
sample_name=${2}

### Initial compute node file system
if [ -d "/tmp/${sample_ID}" ];then
        rm -r "/tmp/${sample_ID}"
        mkdir "/tmp/${sample_ID}"
else
        mkdir "/tmp/${sample_ID}"
fi
if [ -d "/tmp/${sample_ID}/db" ];then
        rm -r "/tmp/${sample_ID}/db"
        mkdir "/tmp/${sample_ID}/db"
else
        mkdir "/tmp/${sample_ID}/db"
fi
if [ -d "/tmp/${sample_ID}/mgf" ];then
        rm -r "/tmp/${sample_ID}/mgf"
        mkdir "/tmp/${sample_ID}/mgf"
else
        mkdir "/tmp/${sample_ID}/mgf"
fi
if [ -d "/tmp/${sample_ID}/result" ];then
        rm -r "/tmp/${sample_ID}/result"
        mkdir "/tmp/${sample_ID}/result"
else
        mkdir "/tmp/${sample_ID}/result"
fi

### Download reference database
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/customprodbj_results/${sample_ID}/protein.pro-ref.fasta /tmp/${sample_ID}/db
### Download mgf file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free/mgf_files/${sample_name}.mgf /tmp/${sample_ID}/mgf
### Download variant peptides
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/rt_prediction/pepquery_input/${sample_ID}_variant_peptides.txt /tmp/${sample_ID}/db

### Run pepquery
java -Xmx60g -jar ~/software/pepquery/pepquery-1.1-jar-with-dependencies.jar \
	-pep /tmp/${sample_ID}/db/${sample_ID}_variant_peptides.txt \
	-db /tmp/${sample_ID}/db/protein.pro-ref.fasta \
	-ms /tmp/${sample_ID}/mgf/${sample_name}.mgf \
	-fixMod 6 \
	-varMod 107 \
	-cpu 8 \
	-minScore 12 \
	-tol 20 \
	-itol 0.05 \
	-n 10000 \
	-um \
	-m 1 \
	-prefix ${sample_ID} \
	-o /tmp/${sample_ID}/result
### Upload results
aws s3 cp /tmp/${sample_ID}/result s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/rt_prediction/pepquery_results/${sample_ID}/ --recursive

### Clear node
rm -r /tmp/${sample_ID}
