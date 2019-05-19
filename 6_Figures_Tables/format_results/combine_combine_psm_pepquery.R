library(tidyverse)

args <- commandArgs(TRUE)
input_file <- args[1]
pepquery <- args[2]
out_put <- args[3]

input_data <- read.delim(input_file)
pepquery_data <- read.delim(pepquery) %>%
        filter(n_db == 0, pvalue <= 0.01, n_ptm == 0, rank == 1)
input_data$index <- NULL
colnames(input_data)[6] <- "y"
colnames(input_data)[3] <- "score"
input_data <- input_data[c(1,2,3,5,6,7,4)]

input_data$pepquery <- ifelse(input_data$peptide %in% pepquery_data$peptide, 1, 0)
write.table(input_data, out_put, row.names = FALSE, quote = FALSE, sep = "\t")

