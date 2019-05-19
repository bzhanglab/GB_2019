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
if [ -d "/tmp/${batch_ID}/somatic_vcf" ];then
        rm -r "/tmp/${batch_ID}/somatic_vcf"
        mkdir "/tmp/${batch_ID}/somatic_vcf"
else
        mkdir "/tmp/${batch_ID}/somatic_vcf"
fi
if [ -d "/tmp/${batch_ID}/germline_vcf" ];then
        rm -r "/tmp/${batch_ID}/germline_vcf"
        mkdir "/tmp/${batch_ID}/germline_vcf"
else
        mkdir "/tmp/${batch_ID}/germline_vcf"
fi
if [ -d "/tmp/${batch_ID}/annovar_results" ];then
        rm -r "/tmp/${batch_ID}/annovar_results"
        mkdir "/tmp/${batch_ID}/annovar_results"
else
        mkdir "/tmp/${batch_ID}/annovar_results"
fi
if [ -d "/tmp/${batch_ID}/annovar_results/somatic" ];then
        rm -r "/tmp/${batch_ID}/annovar_results/somatic"
        mkdir "/tmp/${batch_ID}/annovar_results/somatic"
else
        mkdir "/tmp/${batch_ID}/annovar_results/somatic"
fi
if [ -d "/tmp/${batch_ID}/annovar_results/germline" ];then
        rm -r "/tmp/${batch_ID}/annovar_results/germline"
        mkdir "/tmp/${batch_ID}/annovar_results/germline"
else
        mkdir "/tmp/${batch_ID}/annovar_results/germline"
fi
if [ -d "/tmp/${batch_ID}/customprodbj_results" ];then
        rm -r "/tmp/${batch_ID}/customprodbj_results"
        mkdir "/tmp/${batch_ID}/customprodbj_results"
else
        mkdir "/tmp/${batch_ID}/customprodbj_results"
fi

### Mapping file
table="/home/kail/projects/2019_04_genome_biology/conf/data_table/CPTAC_COprospective_Sample_TMT_mapping_file.txt"
### Input file for customprodbj
input_file="/tmp/${batch_ID}/annovar_results/${batch_ID}_germline_somatic.txt"
echo "sample    somatic germline" >> $input_file
sed 1d $table | while read -r sampleID tmtID
do
	if [ "$tmtID" == "$batch_ID" ]
        then
		echo "${sampleID}	/input/annovar_results/somatic/${sampleID}.hg19_multianno.txt	/input/annovar_results/germline/${sampleID}.hg19_multianno.txt" >> $input_file
		### Download somatic vcf file from S3
		aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/vcf_files/somatic_vcf/${sampleID}.vcf /tmp/${batch_ID}/somatic_vcf/
		### Convert vcf file format
		perl /home/kail/software/annovar/convert2annovar.pl -format vcf4 /tmp/${batch_ID}/somatic_vcf/${sampleID}.vcf > /tmp/${batch_ID}/somatic_vcf/${sampleID}.vcf.avinput
		### Variants annotation by annovar
		perl /home/kail/software/annovar/table_annovar.pl /tmp/${batch_ID}/somatic_vcf/${sampleID}.vcf.avinput /home/kail/projects/2019_04_genome_biology/conf/reference/hg19_mod/ \
		-buildver hg19 \
		-out /tmp/${batch_ID}/annovar_results/somatic/${sampleID} \
		-protocol refGene \
                -operation g \
                -nastring . \
                --thread 8 \
                --maxgenethread 8 \
                -polish \
                -otherinfo
		
		### Download germline vcf file from S3
		aws s3 cp s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/vcf_files/germline_vcf/${sampleID}.vcf /tmp/${batch_ID}/germline_vcf/
		### Convert vcf file format
		perl /home/kail/software/annovar/convert2annovar.pl -format vcf4 /tmp/${batch_ID}/germline_vcf/${sampleID}.vcf > /tmp/${batch_ID}/germline_vcf/${sampleID}.vcf.avinput
		### Variants annotation by annovar
		perl /home/kail/software/annovar/table_annovar.pl /tmp/${batch_ID}/germline_vcf/${sampleID}.vcf.avinput /home/kail/projects/2019_04_genome_biology/conf/reference/hg19_mod/ \
		-buildver hg19 \
		-out /tmp/${batch_ID}/annovar_results/germline/${sampleID} \
		-protocol refGene \
                -operation g \
                -nastring . \
                --thread 8 \
                --maxgenethread 8 \
                -polish \
                -otherinfo
	fi
done
### Upload annovar results to S3 
aws s3 cp /tmp/${batch_ID}/annovar_results/ s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/annovar_results/${batch_ID}/ --recursive
### Build database by customprodbj
docker run -u 510 -v /home/kail/projects/2019_04_genome_biology/conf/reference/hg19_mod/:/database/ -v /tmp/${batch_ID}/:/input/ -v /tmp/${batch_ID}/customprodbj_results/:/output/ zhanglab18/customprodbj:1.1.0 java -jar customprodbj.jar \
	-f /input/annovar_results/${batch_ID}_germline_somatic.txt \
	-d /database/hg19_refGeneMrna.fa \
	-r /database/hg19_refGene.txt \
	-ref /output/protein.pro-ref.fasta \
	-t \
	-o /output/
### Remove redundant protein information in fasta db file to avoid bug in next analysis
java -cp /home/kail/projects/2019_04_genome_biology/conf/scripts/ thanks /tmp/${batch_ID}/customprodbj_results/${batch_ID}_germline_somatic-var.fasta
### Upload new database
aws s3 cp /tmp/${batch_ID}/customprodbj_results/ s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/customprodbj_results/${batch_ID}/ --recursive
### Generate target_decoy database by pga
Rscript pga_generate_decoy_db.R /tmp/${batch_ID}/customprodbj_results/${batch_ID}_germline_somatic-var.fasta.new.fasta /tmp/${batch_ID}/customprodbj_results/${batch_ID}_target_decoy.fasta
### Split target_decoy database to normal and variant for two step search
java -cp /home/kail/projects/2019_04_genome_biology/conf/softwares/mgfIndex-1.0-SNAPSHOT/mgfIndex-1.0-SNAPSHOT.jar splitFasta /tmp/${batch_ID}/customprodbj_results/${batch_ID}_target_decoy.fasta /tmp/${batch_ID}/customprodbj_results/
### Upload the results to S3
aws s3 cp /tmp/${batch_ID}/customprodbj_results/ s3://zhanglab-kail/projects/2019_04_genome_biology/prospective_colon/custom_database/${batch_ID}/ --recursive

### Clear node
rm /tmp/${batch_ID}
