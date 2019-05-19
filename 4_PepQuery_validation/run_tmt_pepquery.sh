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

batch_ID=${1}

### Initial compute node file system
if [ -d "/tmp/${batch_ID}" ];then
        rm -r "/tmp/${batch_ID}"
        mkdir "/tmp/${batch_ID}"
else
        mkdir "/tmp/${batch_ID}"
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
if [ -d "/tmp/${batch_ID}/result" ];then
        rm -r "/tmp/${batch_ID}/result"
        mkdir "/tmp/${batch_ID}/result"
else
        mkdir "/tmp/${batch_ID}/result"
fi
### Download reference database
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/customprodbj_results/${batch_ID}/protein.pro-ref.fasta /tmp/${batch_ID}/db
### Download mgf file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/${batch_ID}.mgf /tmp/${batch_ID}/mgf
### Download variant peptides
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/pepquery_input/${batch_ID}_variant_peptides.txt /tmp/${batch_ID}/db

### Run pepquery
java -Xmx60g -jar ~/software/pepquery/pepquery-1.1-jar-with-dependencies.jar \
	-pep /tmp/${batch_ID}/db/${batch_ID}_variant_peptides.txt \
	-db /tmp/${batch_ID}/db/protein.pro-ref.fasta \
	-ms /tmp/${batch_ID}/mgf/${batch_ID}.mgf \
	-fixMod 6,55,98 \
	-varMod 107 \
	-cpu 8 \
	-minScore 12 \
	-tol 20 \
	-itol 0.05 \
	-n 10000 \
	-um \
	-m 1 \
	-prefix ${batch_ID} \
	-o /tmp/${batch_ID}/result
### Upload results
aws s3 cp /tmp/${batch_ID}/result s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/pepquery_results/${batch_ID}/ --recursive

### Clear node
rm -r /tmp/${batch_ID}
