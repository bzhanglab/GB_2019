#!/bin/sh

#SBATCH --job-name=2nd_step
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

### Download mgf file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/${batch_ID}.mgf /tmp/${batch_ID}/mgf
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/${batch_ID}.mgf.cui /tmp/${batch_ID}/mgf
### Download mgf index file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/index_files/${batch_ID}.mgf_index.txt /tmp/${batch_ID}/mgf
### Download 1st step result
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/${batch_ID}_1st.mzid /tmp/${batch_ID}/1st_step/
### Download database file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/custom_database/${batch_ID}/${batch_ID}_only_normal.fasta /tmp/${batch_ID}/db
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/custom_database/${batch_ID}/${batch_ID}_only_variant.fasta /tmp/${batch_ID}/db

### Convert xml to mzid by mzidentml-lib
java -Xmx60g -jar /home/kail/projects/2019_04_genome_biology/conf/softwares/mzidlib/mzidentml-lib-1.6.12-SNAPSHOT.jar Tandem2mzid /tmp/${batch_ID}/1st_step/${batch_ID}_1st.xml /tmp/${batch_ID}/1st_step/${batch_ID}_1st.mzid -outputFragmentation false -decoyRegex XXX_ -databaseFileFormatID MS:1001348 -massSpecFileFormatID MS:1001062 -idsStartAtZero false -compress false -proteinCodeRegex \S+
aws s3 cp /tmp/${batch_ID}/1st_step/${batch_ID}_1st.mzid s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/1st_step/
### Get 1st spectrum title
Rscript generate_1st_psm_txt.R ${batch_ID} /tmp/${batch_ID}
aws s3 cp /tmp/${batch_ID}/1st_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/1st_step/ --recursive
### Generate 2nd step mgf 
java -Xmx60g -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar getSecondMGF /tmp/${batch_ID}/mgf/${batch_ID}.mgf /tmp/${batch_ID}/2nd_step/1st_step_spectrum_title.txt /tmp/${batch_ID}/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${batch_ID}/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/2nd_step/ --recursive

### Run 2nd step
/home/kail/software/xtandem/tandem.exe /home/kail/projects/2019_04_genome_biology/conf/scripts/xtandem/parameters/${batch_ID}_2nd_par.xml
aws s3 cp /tmp/${batch_ID}/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/xtandem/${batch_ID}/2nd_step/ --recursive

### Clear node
rm -r /tmp/${batch_ID}
