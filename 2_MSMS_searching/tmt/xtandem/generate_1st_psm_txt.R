library(PGA)
library(tidyverse)

args <- commandArgs(TRUE)
batch_ID=args[1]
work_dir=args[2]

parserGear(file=paste(work_dir, "/1st_step/", batch_ID, "_1st.mzid", sep=""),
	db=paste(work_dir, "/db/", batch_ID, "_only_normal.fasta", sep=""),
	decoyPrefix="XXX_",
	prefix=paste(batch_ID, "_1st", sep=""),
	xmx=60,
        thread=8,
        alignment=0,
	outdir=paste(work_dir, "/1st_step/", sep=""))

raw_psm <- read.delim(paste(work_dir, "/1st_step/", batch_ID, "_1st-rawPSMs.txt", sep=""))
colnames(raw_psm)[2] <- "score"
write.table(raw_psm, paste(work_dir, "/1st_step/", batch_ID, "_1st-rawPSMs.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")

calculateFDR(psmfile=paste(work_dir, "/1st_step/", batch_ID, "_1st-rawPSMs.txt", sep=""),
	db=paste(work_dir, "/db/", batch_ID, "_only_normal.fasta", sep=""),
	fdr=0.01,
        decoyPrefix="XXX_",
        better_score_lower=FALSE,
        remap=FALSE,
        score_t = 0,
        protein_inference=FALSE,
	out_dir=paste(work_dir, "/1st_step/", sep=""),
	xmx=60)

pga_global_data <- read.delim(paste(work_dir, "/1st_step/pga-peptideSummary.txt", sep=""))
mgf_index_data <- read.delim(paste(work_dir, "/mgf/", batch_ID, ".mgf_index.txt", sep=""))
output_file <- paste(work_dir, "/2nd_step/1st_step_spectrum_title.txt", sep="")
filter_psm_data <- pga_global_data %>%
        filter(Qvalue <= 0.01) %>%
	filter(str_length(peptide) > 6, str_length(peptide) < 46) %>%
	filter(!str_detect(mods, "42.0105"), !str_detect(mods, "18.0105"), !str_detect(mods, "17.02")) %>%
        mutate(index = as.numeric(index) + 1)
pga_global_data_with_fraction <- left_join(filter_psm_data, mgf_index_data, by = "index")
write.table(pga_global_data_with_fraction$Title, output_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
