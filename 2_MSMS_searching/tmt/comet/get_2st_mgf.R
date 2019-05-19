library(tidyverse)

args <- commandArgs(TRUE)
pga_global_file <- args[1]
mgf_index_file <- args[2]
output_file <- args[3]

pga_global_data <- read.delim(pga_global_file)
mgf_index_data <- read.delim(mgf_index_file)
filter_psm_data <- pga_global_data %>%
        filter(Qvalue <= 0.01) %>%
	filter(str_length(peptide) > 6, str_length(peptide) < 46) %>%
        mutate(index = as.numeric(index))
pga_global_data_with_fraction <- left_join(filter_psm_data, mgf_index_data, by = "index")
write.table(pga_global_data_with_fraction$Title, output_file, row.names = FALSE, quote = FALSE, sep = "\t", col.names = FALSE)
