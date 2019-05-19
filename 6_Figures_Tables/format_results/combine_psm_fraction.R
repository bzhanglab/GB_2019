library(tidyverse)

args <- commandArgs(TRUE)
input_file <- args[1]
mgf_index_file <- args[2]
out_put_file <- args[3]

combined_data <- read.delim(input_file)
mgf_index_data <- read.delim(mgf_index_file)
filter_var_psm_data <- combined_data%>%
  filter(Qvalue <= 0.01, !str_detect(protein, "NM"), !str_detect(protein, "cont_")) %>%
        filter(str_length(peptide) > 6, str_length(peptide) < 46)
process_var_data <- filter_var_psm_data %>%
  select(index, peptide, mods, evalue)
process_var_data$file <- str_sub(process_var_data$index, 1, 44)
process_var_data$index <- str_sub(process_var_data$index, 46)
process_var_data <- process_var_data %>%
  mutate(index = as.numeric(index) + 1)
var_data_with_fraction <- left_join(process_var_data, mgf_index_data, by = c("file","index"))
var_data_with_fraction$fraction <- str_sub(var_data_with_fraction$Title, 1, 46)
var_data_with_fraction$rt <- var_data_with_fraction$rt/60
colnames(var_data_with_fraction)[5] <- "sample"
write.table(var_data_with_fraction, out_put_file, row.names = FALSE, quote = FALSE, sep = "\t")
