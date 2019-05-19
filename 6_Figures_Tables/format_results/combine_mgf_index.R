library(tidyverse)

args <- commandArgs(TRUE)
index_folder <- args[1]

index_files <- dir(index_folder, pattern="Proteome")
all_index_data <- data_frame()
for (file in index_files) {
	tmt_name <- str_replace_all(file, ".mgf_index.txt", "")
	one_data <- read.delim(paste(index_folder, file, sep=""))
	one_data$file <- tmt_name
	all_index_data <- bind_rows(all_index_data, one_data)
}
write.table(all_index_data, paste(index_folder, "all_index_information.txt", sep=""), row.names = FALSE, quote = FALSE, sep = "\t")

