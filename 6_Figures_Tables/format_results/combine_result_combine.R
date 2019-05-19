library(tidyverse)

args <- commandArgs(TRUE)
global <- args[1]
separate <- args[2]
two_step <- args[3]
pepquery <- args[4]
out_put <- args[5]

global_data <- read.delim(global)
colnames(global_data)[7] <- "global_rt"

separate_data <- read.delim(separate)
colnames(separate_data)[7] <- "separate_rt"
part_separate_data <- separate_data %>%
        select(peptide, separate_rt)
two_step_data_or <- read.delim(two_step)
two_step_data <- two_step_data_or[order(two_step_data_or$peptide, -(two_step_data_or$evalue)), ]
two_step_data <- two_step_data[ !duplicated(two_step_data$peptide), ]
colnames(two_step_data)[6] <- "two_step_rt"
part_two_step_data <- two_step_data %>%
        select(peptide, two_step_rt)

pepquery_data <- read.delim(pepquery) %>%
        filter(n_db == 0, pvalue <= 0.01, n_ptm == 0, rank == 1)

combine_data <- global_data
combine_data$global <- ifelse(combine_data$peptide %in% global_data$peptide, 1, 0)
combine_data$separate <- ifelse(combine_data$peptide %in% separate_data$peptide, 1, 0)
combine_data$two_step <- ifelse(combine_data$peptide %in% two_step_data$peptide, 1, 0)

sep_diff <- setdiff(separate_data$peptide, combine_data$peptide)
if (length(sep_diff) != 0){
        sep_only <- separate_data %>%
                filter(peptide %in% sep_diff)
	combine_data <- left_join(combine_data, part_separate_data, by="peptide")
        combine_data <- bind_rows(combine_data, sep_only)
        combine_data$separate <- ifelse(combine_data$peptide %in% separate_data$peptide, 1, 0)
        combine_data$two_step <- ifelse(combine_data$peptide %in% two_step_data$peptide, 1, 0)
        combine_data$global <- ifelse(combine_data$peptide %in% global_data$peptide, 1, 0)
} else{
        combine_data <- left_join(combine_data, part_separate_data, by="peptide")
}
two_diff <- setdiff(two_step_data$peptide, combine_data$peptide)
if (length(two_diff) != 0){
	two_only <- two_step_data %>%
                filter(peptide %in% two_diff)
        combine_data <- left_join(combine_data, part_two_step_data, by="peptide")
        combine_data <- bind_rows(combine_data, two_only)
        combine_data$separate <- ifelse(combine_data$peptide %in% separate_data$peptide, 1, 0)
        combine_data$two_step <- ifelse(combine_data$peptide %in% two_step_data$peptide, 1, 0)
        combine_data$global <- ifelse(combine_data$peptide %in% global_data$peptide, 1, 0)
} else{
        combine_data <- left_join(combine_data, part_two_step_data, by="peptide")
}
combine_data$pepquery <- ifelse(combine_data$peptide %in% pepquery_data$peptide, 1, 0)

process_combine_data <- combine_data[order(combine_data$peptide, -(combine_data$evalue)), ]
process_combine_data <- process_combine_data[ !duplicated(process_combine_data$peptide), ]

write.table(combine_data, out_put, row.names = FALSE, quote = FALSE, sep = "\t")

