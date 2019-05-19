library(PGA)

args <- commandArgs(TRUE)
input=args[1]
output=args[2]

buildTargetDecoyDB(db=input,
	decoyPrefix="XXX_",
	cont_file="/home/kail/projects/2019_04_genome_biology/conf/database/contaminants.fasta",
	output=output)
