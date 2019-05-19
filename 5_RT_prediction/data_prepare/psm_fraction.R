library(tidyverse)

args <- commandArgs(TRUE)
pga_global_file <- args[1]
mgf_index_file <- args[2]
output_file <- args[3]

pga_global_data <- read.delim(pga_global_file)
mgf_index_data <- read.delim(mgf_index_file)
filter_psm_data <- pga_global_data %>%
  filter(Qvalue <= 0.01, str_detect(protein, "NM")) 
#filter_psm_data <- filter_psm_data[order(filter_psm_data$peptide, (filter_psm_data$specEValue)), ]
#filter_psm_data <- filter_psm_data[ !duplicated(filter_psm_data$peptide), ]

filter_var_psm_data <- pga_global_data %>%
  filter(Qvalue <= 0.01, !str_detect(protein, "NM"), !str_detect(protein, "cont_"))
#filter_var_psm_data <- filter_var_psm_data[order(filter_var_psm_data$peptide, (filter_var_psm_data$specEValue)), ]
#filter_var_psm_data <- filter_var_psm_data[ !duplicated(filter_var_psm_data$peptide), ]

process_data <- filter_psm_data %>%
  select(index, peptide, mods, evalue, rt) %>%
  mutate(index = as.numeric(index) + 1)
process_data$rt <- process_data$rt/60

process_var_data <- filter_var_psm_data %>%
  select(index, peptide, mods, evalue, rt) %>%
  mutate(index = as.numeric(index) + 1)
process_var_data$rt <- process_var_data$rt/60

pga_global_data_with_fraction <- left_join(process_data, mgf_index_data, by = "index")
pga_global_var_data_with_fraction <- left_join(process_var_data, mgf_index_data, by = "index")
write.table(pga_global_data_with_fraction, paste(output_file, "_normal_psm.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")
write.table(pga_global_var_data_with_fraction, paste(output_file, "_variant_psm.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")
