library(PGA)
library(tidyverse)

args <- commandArgs(TRUE)
input <- args[1]
output_path <- args[2]
decoy_tag <- args[3]
database <- args[4]

raw_psm <- read.delim(input)
write.table(raw_psm, input, row.names = FALSE, quote = FALSE, sep = "\t")

calculateFDR(psmfile=input,
        db=database,
        fdr=0.01,
        decoyPrefix=decoy_tag,
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=TRUE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output_path, "/global_fdr", sep=""),
        xmx=8)
calculateFDR(psmfile=input,
        db=database,
        fdr=0.01,
        decoyPrefix=decoy_tag,
        novelPrefix="VAR",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=TRUE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output_path, "/separate_fdr", sep=""),
        xmx=8)
