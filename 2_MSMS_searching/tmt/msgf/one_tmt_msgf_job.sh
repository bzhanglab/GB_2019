#!/bin/sh

#SBATCH --job-name=colon_msgf
#SBATCH --output=%A_%a.out
#SBATCH --error=%A_%a.err
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=50000M
#SBATCH --nodes=1
#SBATCH --ntasks=1

set -e -x

batch_ID=${1}

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
if [ -d "/tmp/${batch_ID}/results/target_decoy" ];then
        rm -r "/tmp/${batch_ID}/results/target_decoy"
        mkdir "/tmp/${batch_ID}/results/target_decoy"
else
        mkdir "/tmp/${batch_ID}/results/target_decoy"
fi
if [ -d "/tmp/${batch_ID}/results/target_decoy/mzid" ];then
        rm -r "/tmp/${batch_ID}/results/target_decoy/mzid"
        mkdir "/tmp/${batch_ID}/results/target_decoy/mzid"
else
        mkdir "/tmp/${batch_ID}/results/target_decoy/mzid"
fi

if [ -d "/tmp/${batch_ID}/results/target_decoy/tsv" ];then
        rm -r "/tmp/${batch_ID}/results/target_decoy/tsv"
        mkdir "/tmp/${batch_ID}/results/target_decoy/tsv"
else
        mkdir "/tmp/${batch_ID}/results/target_decoy/tsv"
fi
if [ -d "/tmp/${batch_ID}/results/normal" ];then
        rm -r "/tmp/${batch_ID}/results/normal"
        mkdir "/tmp/${batch_ID}/results/normal"
else
        mkdir "/tmp/${batch_ID}/results/normal"
fi
if [ -d "/tmp/${batch_ID}/results/normal/mzid" ];then
        rm -r "/tmp/${batch_ID}/results/normal/mzid"
        mkdir "/tmp/${batch_ID}/results/normal/mzid"
else
        mkdir "/tmp/${batch_ID}/results/normal/mzid"
fi

if [ -d "/tmp/${batch_ID}/results/normal/tsv" ];then
        rm -r "/tmp/${batch_ID}/results/normal/tsv"
        mkdir "/tmp/${batch_ID}/results/normal/tsv"
else
        mkdir "/tmp/${batch_ID}/results/normal/tsv"
fi

aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/custom_database/${batch_ID}/ /tmp/${batch_ID}/db --recursive
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/mgf_files/${batch_ID}.mgf /tmp/${batch_ID}/mgf
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/msgf+/Mod_files/Mods.txt /tmp/${batch_ID}/db

#java -Xmx48g -jar ~/software/MSGFPlus_v20181015/MSGFPlus.jar \
#	-s /tmp/${batch_ID}/mgf/${batch_ID}.mgf \
#	-d /tmp/${batch_ID}/db/${batch_ID}_with_cont.fasta \
#	-mod /tmp/${batch_ID}/db/Mods.txt \
#	-o /tmp/${batch_ID}/results/normal/mzid/${batch_ID}.mzid \
#	-m 3 \
#	-maxLength 45 \
#	-n 1 \
#	-minlength 7 \
#	-addFeatures 1 \
#	-t 20.0ppm \
#	-minCharge 2 \
#	-maxCharge 3 \
#	-tda 1 \
#	-protocol 4 \
#	-ntt 2 \
#	-inst 1 \
#	-ti 0,0 \
#	-e 1

#java -Xmx48g -jar ~/software/Deeplfq/DeepLFQ-1.0.0-jar-with-dependencies.jar \
#	-d /tmp/${batch_ID}/db/${batch_ID}_with_cont.fasta \
#	-i /tmp/${batch_ID}/results/normal/mzid \
#	-o /tmp/${batch_ID}/results/normal/tsv

#java -Xmx48g -jar ~/software/MSGFPlus_v20181015/MSGFPlus.jar \
#        -s /tmp/${batch_ID}/mgf/${batch_ID}.mgf \
#        -d /tmp/${batch_ID}/db/${batch_ID}_target_decoy.fasta \
#        -mod /tmp/${batch_ID}/db/Mods.txt \
#        -o /tmp/${batch_ID}/results/target_decoy/mzid/${batch_ID}.mzid \
#        -m 3 \
#        -maxLength 45 \
#        -n 1 \
#        -minlength 7 \
#        -addFeatures 1 \
#        -t 20.0ppm \
#        -minCharge 2 \
#        -maxCharge 3 \
#        -tda 0 \
#        -protocol 4 \
#        -ntt 2 \
#        -inst 1 \
#        -ti 0,0 \
#        -e 1

java -Xmx48g -jar ~/software/MSGFPlus_v20181015/MSGFPlus.jar \
        -s /tmp/${batch_ID}/mgf/${batch_ID}.mgf \
        -d /tmp/${batch_ID}/db/${batch_ID}_only_normal.fasta \
        -mod /tmp/${batch_ID}/db/Mods.txt \
        -o /tmp/${batch_ID}/results/${batch_ID}_1st.mzid \
        -m 3 \
        -maxLength 45 \
        -n 1 \
        -minlength 7 \
        -addFeatures 1 \
        -t 20.0ppm \
        -minCharge 2 \
        -maxCharge 3 \
        -tda 0 \
        -protocol 4 \
        -ntt 2 \
        -inst 1 \
        -ti 0,0 \
        -e 1

#java -Xmx48g -jar ~/software/Deeplfq/DeepLFQ-1.0.0-jar-with-dependencies.jar \
#        -d /tmp/${batch_ID}/db/${batch_ID}_target_decoy.fasta \
#        -i /tmp/${batch_ID}/results/target_decoy/mzid \
#        -o /tmp/${batch_ID}/results/target_decoy/tsv

aws s3 cp /tmp/${batch_ID}/results/${batch_ID}_1st.mzid s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/proteome_results/${batch_ID}/
rm -r "/tmp/${batch_ID}"
