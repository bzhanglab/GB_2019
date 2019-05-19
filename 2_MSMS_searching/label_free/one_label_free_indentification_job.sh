#!/bin/sh

#SBATCH --job-name=label_free
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50000M
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -e -x

batch_ID=${1}
sample_ID=${2}
sample_name=${3}

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
if [ -d "/tmp/${sample_ID}/results" ];then
        rm -r "/tmp/${sample_ID}/results"
        mkdir "/tmp/${sample_ID}/results"
else
        mkdir "/tmp/${sample_ID}/results"
fi
if [ -d "/tmp/${sample_ID}/results/msgf" ];then
        rm -r "/tmp/${sample_ID}/results/msgf"
        mkdir "/tmp/${sample_ID}/results/msgf"
else
        mkdir "/tmp/${sample_ID}/results/msgf"
fi
if [ -d "/tmp/${sample_ID}/results/comet" ];then
        rm -r "/tmp/${sample_ID}/results/comet"
        mkdir "/tmp/${sample_ID}/results/comet"
else
        mkdir "/tmp/${sample_ID}/results/comet"
fi
if [ -d "/tmp/${sample_ID}/results/xtandem" ];then
        rm -r "/tmp/${sample_ID}/results/xtandem"
        mkdir "/tmp/${sample_ID}/results/xtandem"
else
        mkdir "/tmp/${sample_ID}/results/xtandem"
fi
### Download mgf file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free/mgf_files/${sample_name}.mgf /tmp/${sample_ID}/mgf/
mv /tmp/${sample_ID}/mgf/${sample_name}.mgf /tmp/${sample_ID}/mgf/${sample_ID}.mgf
### Download mgf index file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free/mgf_files/index_files/${sample_name}.mgf_index.txt /tmp/${sample_ID}/mgf/
### Download database file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/customprodbj_results/${batch_ID}/${sample_ID}/ /tmp/${sample_ID}/db/ --recursive
### Download msgf modification file
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/msgf+/Mod_files/label_free_Mods.txt /tmp/${sample_ID}/db

### Run msgf+
./run_label_free_msgf.sh /tmp/${sample_ID}/mgf/${sample_ID}.mgf /tmp/${sample_ID}/db/${sample_ID}_target_decoy.fasta /tmp/${sample_ID}/results/msgf/${sample_ID}.mzid /tmp/${sample_ID}/db/label_free_Mods.txt
### Run msgf+ 2 step 1st step
./run_label_free_msgf.sh /tmp/${sample_ID}/mgf/${sample_ID}.mgf /tmp/${sample_ID}/db/${sample_ID}_only_normal.fasta /tmp/${sample_ID}/results/msgf/${sample_ID}_1st.mzid /tmp/${sample_ID}/db/label_free_Mods.txt
### Upload results to S3
aws s3 cp /tmp/${sample_ID}/results/msgf s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/msgf/${sample_ID}/ --recursive

### Caculate global and separate fdr
mkdir /tmp/${sample_ID}/results/msgf/global_fdr
mkdir /tmp/${sample_ID}/results/msgf/separate_fdr
Rscript caculate_fdr.R /tmp/${sample_ID}/results/msgf/${sample_ID}.mzid /tmp/${sample_ID}/db/${sample_ID}_target_decoy.fasta /tmp/${sample_ID}/results/msgf/ ${sample_ID}
aws s3 cp /tmp/${sample_ID}/results/msgf s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/msgf/${sample_ID}/ --recursive
### Caculate msgf 1st step fdr
mkdir /tmp/${sample_ID}/results/msgf/1st_step
mv /tmp/${sample_ID}/results/msgf/${sample_ID}_1st.mzid /tmp/${sample_ID}/results/msgf/1st_step/
mkdir /tmp/${sample_ID}/results/msgf/2nd_step
Rscript generate_1st_psm_txt.R ${sample_ID} /tmp/${sample_ID}/results/msgf/ /tmp/${sample_ID}/ $sample_name
aws s3 cp /tmp/${sample_ID}/results/msgf/1st_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/msgf/${sample_ID}/1st_step/ --recursive
### Generate msgf 2nd mgf
java -Xmx60g -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar getSecondMGF /tmp/${sample_ID}/mgf/${sample_ID}.mgf /tmp/${sample_ID}/results/msgf/2nd_step/1st_step_spectrum_title.txt /tmp/${sample_ID}/results/msgf/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${sample_ID}/results/msgf/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/msgf/${sample_ID}/2nd_step/ --recursive
### Run msgf second step
./run_msgf.sh /tmp/${sample_ID}/results/msgf/2nd_step/2nd_step.mgf /tmp/${sample_ID}/db/${sample_ID}_only_variant.fasta /tmp/${sample_ID}/results/msgf/2nd_step/${sample_ID}_2nd.mzid /tmp/${sample_ID}/db/label_free_Mods.txt
aws s3 cp /tmp/${sample_ID}/results/msgf/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/msgf/${sample_ID}/2nd_step/ --recursive
### Caculate 2 step 2nd step fdr by pga
Rscript generate_2nd_psm_txt.R ${sample_ID} /tmp/${sample_ID}/results/msgf/ /tmp/${sample_ID}/
aws s3 cp /tmp/${sample_ID}/results/msgf/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/msgf/${sample_ID}/2nd_step/ --recursive
### Run xtandem
mkdir /tmp/${sample_ID}/results/xtandem/1st_step
mkdir /tmp/${sample_ID}/results/xtandem/2nd_step
./run_xtandem.sh ${sample_ID}
aws s3 cp /tmp/${sample_ID}/results/xtandem s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/ --recursive
### Convert xml to mzid by mzidentml-lib
java -Xmx60g -jar /home/kail/projects/2019_04_genome_biology/conf/softwares/mzidlib/mzidentml-lib-1.6.12-SNAPSHOT.jar Tandem2mzid /tmp/${sample_ID}/results/xtandem/${sample_ID}_decoy.xml /tmp/${sample_ID}/results/xtandem/${sample_ID}.mzid -outputFragmentation false -decoyRegex XXX_ -databaseFileFormatID MS:1001348 -massSpecFileFormatID MS:1001062 -idsStartAtZero false -compress false -proteinCodeRegex "\\S+"
java -Xmx60g -jar /home/kail/projects/2019_04_genome_biology/conf/softwares/mzidlib/mzidentml-lib-1.6.12-SNAPSHOT.jar Tandem2mzid /tmp/${sample_ID}/results/xtandem/1st_step/${sample_ID}_1st.xml /tmp/${sample_ID}/results/xtandem/1st_step/${sample_ID}_1st.mzid -outputFragmentation false -decoyRegex XXX_ -databaseFileFormatID MS:1001348 -massSpecFileFormatID MS:1001062 -idsStartAtZero false -compress false -proteinCodeRegex "\\S+"
aws s3 cp /tmp/${sample_ID}/results/xtandem s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/ --recursive
### Caculate xtandem results fdr
mkdir /tmp/${sample_ID}/results/xtandem/global_fdr
mkdir /tmp/${sample_ID}/results/xtandem/separate_fdr
Rscript caculate_fdr.R /tmp/${sample_ID}/results/xtandem/${sample_ID}.mzid /tmp/${sample_ID}/db/${sample_ID}_target_decoy.fasta /tmp/${sample_ID}/results/xtandem/ ${sample_ID}
aws s3 cp /tmp/${sample_ID}/results/xtandem s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/ --recursive
### Caculate Xtandem 1st step fdr
Rscript generate_1st_psm_txt.R ${sample_ID} /tmp/${sample_ID}/results/xtandem/ /tmp/${sample_ID}/ $sample_name
aws s3 cp /tmp/${sample_ID}/results/xtandem/1st_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/1st_step/ --recursive
### Get Xtandem 2nd mgf
java -Xmx60g -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar getSecondMGF /tmp/${sample_ID}/mgf/${sample_ID}.mgf /tmp/${sample_ID}/results/xtandem/2nd_step/1st_step_spectrum_title.txt /tmp/${sample_ID}/results/xtandem/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${sample_ID}/results/xtandem/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/2nd_step/ --recursive
### Run Xtandem second step
/home/kail/software/xtandem/tandem.exe /home/kail/projects/2019_04_genome_biology/conf/scripts/label_free/xtandem_parameters/${sample_ID}_2nd_par.xml
aws s3 cp /tmp/${sample_ID}/results/xtandem/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/2nd_step/ --recursive
java -Xmx60g -jar /home/kail/projects/2019_04_genome_biology/conf/softwares/mzidlib/mzidentml-lib-1.6.12-SNAPSHOT.jar Tandem2mzid /tmp/${sample_ID}/results/xtandem/2nd_step/${sample_ID}_2nd.xml /tmp/${sample_ID}/results/xtandem/2nd_step/${sample_ID}_2nd.mzid -outputFragmentation false -decoyRegex XXX_ -databaseFileFormatID MS:1001348 -massSpecFileFormatID MS:1001062 -idsStartAtZero false -compress false -proteinCodeRegex "\\S+"
aws s3 cp /tmp/${sample_ID}/results/xtandem/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/2nd_step/ --recursive
Rscript generate_2nd_psm_txt.R ${sample_ID} /tmp/${sample_ID}/results/xtandem/ /tmp/${sample_ID}/
aws s3 cp /tmp/${sample_ID}/results/xtandem/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/xtandem/${sample_ID}/2nd_step/ --recursive
### Comet search
mkdir /tmp/${sample_ID}/results/comet/1st_step
mkdir /tmp/${sample_ID}/results/comet/2nd_step
mkdir /tmp/${sample_ID}/results/comet/global_fdr
mkdir /tmp/${sample_ID}/results/comet/separate_fdr
./run_comet.sh ${sample_ID} /tmp/${sample_ID}/results/comet/${sample_ID} /tmp/${sample_ID}/mgf/${sample_ID}.mgf /tmp/${sample_ID}/db /tmp/${sample_ID}/results/comet
aws s3 cp /tmp/${sample_ID}/results/comet s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/comet/${sample_ID}/ --recursive
### Caculate comet result fdr
Rscript pga_fdr_comet_1st.R ${sample_ID} /tmp/${sample_ID}/results/comet/ /tmp/${sample_ID}/ $sample_name
aws s3 cp /tmp/${sample_ID}/results/comet/1st_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/comet/${sample_ID}/1st_step/ --recursive
### Generate comet 2nd mgf
java -Xmx60g -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar getSecondMGF /tmp/${sample_ID}/mgf/${sample_ID}.mgf /tmp/${sample_ID}/results/comet/2nd_step/1st_step_spectrum_title.txt /tmp/${sample_ID}/results/comet/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${sample_ID}/results/comet/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/comet/${sample_ID}/2nd_step/ --recursive
### Run comet second step
/home/kail/software/comet/comet.2018014.linux.exe -P/home/kail/projects/2019_04_genome_biology/conf/scripts/label_free/comet_parameters/${sample_ID}_2nd.params -N/tmp/${sample_ID}/results/comet/2nd_step/${sample_ID}_2nd /tmp/${sample_ID}/results/comet/2nd_step/2nd_step.mgf
aws s3 cp /tmp/${sample_ID}/results/comet/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/comet/${sample_ID}/2nd_step/ --recursive
./format_comet_results.sh /tmp/${sample_ID}/results/comet/2nd_step/${sample_ID}_2nd.txt
### Caculate comet 2nd step fdr
Rscript pga_fdr_2nd_comet.R /tmp/${sample_ID}/results/comet/2nd_step/${sample_ID}_2nd.txt_for_pga.txt /tmp/${sample_ID}/results/comet/2nd_step /tmp/${sample_ID}/db/${sample_ID}_only_variant.fasta
aws s3 cp /tmp/${sample_ID}/results/comet/2nd_step s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/label_free_new_mod/results/comet/${sample_ID}/2nd_step/ --recursive

### Clear node
rm -r /tmp/${sample_ID} 
