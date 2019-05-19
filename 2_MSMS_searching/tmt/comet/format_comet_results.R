library(tidyverse)

####### Comet result has a blank in last column; Add one more column name manually.
args <- commandArgs(TRUE)
input <- args[1]
output <- args[2]

comet_result <- read.delim(input) %>%
        select(scan, "xcorr", charge, exp_neutral_mass, plain_peptide, protein, modifications)
colnames(comet_result) <- c("index", "score", "charge", "mass", "peptide", "protein", "mods")
comet_result$protein <- str_replace_all(comet_result$protein, ",", ";")
comet_result$mods <- str_replace_all(comet_result$mods, ",", ";")

write_tsv(comet_result, output)
