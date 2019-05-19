library(PGA)

args <- commandArgs(TRUE)
input_raw <- args[1]
database <- args[2]
output <- args[3]

calculateFDR(psmfile=input_raw,
        db=database,
        fdr=0.01,
        decoyPrefix="XXX_",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=FALSE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output, "global_fdr", sep=""),
        xmx=20)
calculateFDR(psmfile=input_raw,
        db=database,
        fdr=0.01,
        decoyPrefix="XXX_",
        novelPrefix="VAR",
        better_score_lower=FALSE,
        remap=FALSE,
        peptide_level=FALSE,
        score_t = 0,
        protein_inference=FALSE,
        out_dir=paste(output, "separate_fdr", sep=""),
        xmx=20)
