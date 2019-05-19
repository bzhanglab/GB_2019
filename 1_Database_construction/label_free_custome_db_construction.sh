#!/bin/sh

#SBATCH --job-name=Build_db
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
if [ -d "/tmp/${batch_ID}/annovar_results" ];then
        rm -r "/tmp/${batch_ID}/annovar_results"
        mkdir "/tmp/${batch_ID}/annovar_results"
else
        mkdir "/tmp/${batch_ID}/annovar_results"
fi
if [ -d "/tmp/${batch_ID}/customprodbj_results" ];then
        rm -r "/tmp/${batch_ID}/customprodbj_results"
        mkdir "/tmp/${batch_ID}/customprodbj_results"
else
	mkdir "/tmp/${batch_ID}/customprodbj_results"
fi
### Download annovar results from S3
aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/annovar_results/${batch_ID}/ /tmp/${batch_ID}/annovar_results/ --recursive
### Mapping file
table="/home/kail/projects/2019_04_genome_biology/conf/data_table/CPTAC_COprospective_TMT_label_free_link.txt"
sed 1d $table | while read -r sampleID tmtID mgfName
do
        if [ "$tmtID" == "$batch_ID" ]
        then
		### Input file for customprodbj
		input_file="/tmp/${batch_ID}/annovar_results/${sampleID}_germline_somatic.txt"
		echo "sample    somatic germline" > $input_file
		echo "${sampleID}       /data/annovar_results/somatic/${sampleID}.hg19_multianno.txt   /data/annovar_results/germline/${sampleID}.hg19_multianno.txt" >> $input_file
		mkdir /tmp/${batch_ID}/customprodbj_results/${sampleID}
		### Build database by customprodbj
		docker run -u 510 -v /home/kail/projects/2019_04_genome_biology/conf/reference/hg19_mod/:/database/ -v /tmp/${batch_ID}/:/data/ zhanglab18/customprodbj:1.1.0 java -jar customprodbj.jar \
		-f /data/annovar_results/${sampleID}_germline_somatic.txt \
		-d /database/hg19_refGeneMrna.fa \
        	-r /database/hg19_refGene.txt \
		-f /data/customprodbj_results/${sampleID}/protein.pro-ref.fasta \
		-t \
		-o /data/customprodbj_results/${sampleID}
		
		### Remove redundant protein information in fasta db file to avoid bug in next analysis
		java -cp /home/kail/projects/2019_04_genome_biology/conf/scripts/ thanks /tmp/${batch_ID}/customprodbj_results/${sampleID}/${sampleID}_germline_somatic-var.fasta
		### Upload new database
		aws s3 cp /tmp/${batch_ID}/customprodbj_results/${sampleID} s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/customprodbj_results/${batch_ID}/${sampleID}/ --recursive
		### Generate target_decoy database by pga
		Rscript pga_generate_decoy_db.R /tmp/${batch_ID}/customprodbj_results/${sampleID}/${sampleID}_germline_somatic-var.fasta.new.fasta /tmp/${batch_ID}/customprodbj_results/${sampleID}/${sampleID}_target_decoy.fasta
		### Split target_decoy database to normal and variant for two step search
		java -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar splitFasta /tmp/${batch_ID}/customprodbj_results/${sampleID}/${sampleID}_target_decoy.fasta /tmp/${batch_ID}/customprodbj_results/${sampleID}/${sampleID}
		### Upload the results to S3
		aws s3 cp /tmp/${batch_ID}/customprodbj_results/${sampleID} s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/customprodbj_results/${batch_ID}/${sampleID}/ --recursive
	fi
done

### Clear node
rm -r /tmp/${batch_ID}
