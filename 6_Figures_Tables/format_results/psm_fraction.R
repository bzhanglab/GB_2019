library(tidyverse)

args <- commandArgs(TRUE)
input_file <- args[1]
mgf_index_file <- args[2]
out_put_file <- args[3]
sample_ID <- args[4]

input_data <- read.delim(input_file)
mgf_index_data <- read.delim(mgf_index_file)
filter_var_psm_data <- input_data %>%
  filter(Qvalue <= 0.01, !str_detect(protein, "NM"), !str_detect(protein, "cont_"))
process_var_data <- filter_var_psm_data %>%
  select(index, peptide, mods, evalue) %>%
  mutate(index = as.numeric(index) + 1)
var_data_with_fraction <- left_join(process_var_data, mgf_index_data, by = "index")
var_data_with_fraction$fraction <- str_sub(var_data_with_fraction$Title, 1, 46)
var_data_with_fraction$rt <- var_data_with_fraction$rt/60
var_data_with_fraction$sample <- sample_ID

write.table(var_data_with_fraction, out_put_file, row.names = FALSE, quote = FALSE, sep = "\t")

