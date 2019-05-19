library(tidyverse)

args <- commandArgs(TRUE)
input_file <- args[1]
pepquery <- args[2]
out_put <- args[3]

input_data <- read.delim(input_file)
#input_data <- input_data_or[order(input_data_or$peptide, -(input_data_or$evalue)), ] # This is for two step combined process
#input_data <- input_data[ !duplicated(input_data$peptide), ] # This is for two step combined process
pepquery_data <- read.delim(pepquery) %>%
        filter(n_db == 0, pvalue <= 0.01, n_ptm == 0, rank == 1)
input_data$index <- NULL
colnames(input_data)[5] <- "y"

input_data$pepquery <- ifelse(input_data$peptide %in% pepquery_data$peptide, 1, 0)
write.table(input_data, out_put, row.names = FALSE, quote = FALSE, sep = "\t")

