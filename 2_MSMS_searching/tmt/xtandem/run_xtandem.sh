#!/bin/sh

#SBATCH --job-name=colon_xtandem
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50000M
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
if [ -d "/tmp/${batch_ID}/results" ];then
        rm -r "/tmp/${batch_ID}/results"
        mkdir "/tmp/${batch_ID}/results"
else
        mkdir "/tmp/${batch_ID}/results"
fi

### Download database file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/custom_database/${batch_ID}/ /tmp/${batch_ID}/db --recursive
### Download mgf file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/${batch_ID}.mgf /tmp/${batch_ID}/mgf

### Run xtandem
/home/kail/software/xtandem/tandem.exe /home/kail/projects/2019_04_genome_biology/conf/scripts/xtandem/parameters/${batch_ID}_decoy_par.xml
aws s3 cp /tmp/${batch_ID}/results/${batch_ID}_decoy.xml s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/
### Run 1st step xtandem
/home/kail/software/xtandem/tandem.exe /home/kail/projects/2019_04_genome_biology/conf/scripts/xtandem/parameters/${batch_ID}_1st_par.xml
aws s3 cp /tmp/${batch_ID}/results/${batch_ID}_1st.xml s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/

### Clear node
rm -r /tmp/${batch_ID}
