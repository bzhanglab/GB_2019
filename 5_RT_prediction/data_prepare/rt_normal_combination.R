library(tidyverse)

args <- commandArgs(TRUE)
msgf_file <- args[1]
xtandem_file <- args[2]
comet_file <- args[3]
out_file <- args[4]

msgf <- read.delim(msgf_file)
xtandem <- read.delim(xtandem_file)
comet <- read.delim(comet_file)

#### msgf and xtandem overlap
msgf_xtandem_inner_res <- inner_join(msgf, xtandem, by="x")
msgf_xtandem_inner_res$error <- abs(msgf_xtandem_inner_res$y.x - msgf_xtandem_inner_res$y.y)
msgf_xtandem_inner_res <- msgf_xtandem_inner_res %>%
  filter(error <= 5)
msgf_xtandem_inner_res$y <- (msgf_xtandem_inner_res$y.x + msgf_xtandem_inner_res$y.y)/2
msgf_xtandem_inner_res <- msgf_xtandem_inner_res %>%
  select(x, y)
#### comet and xtandem overlap
comet_xtandem_inner_res <- inner_join(comet, xtandem, by="x")
comet_xtandem_inner_res$error <- abs(comet_xtandem_inner_res$y.x - comet_xtandem_inner_res$y.y)
comet_xtandem_inner_res <- comet_xtandem_inner_res %>%
  filter(error <= 5)
comet_xtandem_inner_res$y <- (comet_xtandem_inner_res$y.x + comet_xtandem_inner_res$y.y)/2
comet_xtandem_inner_res <- comet_xtandem_inner_res %>%
  select(x, y)
#### msgf and comet overlap
msgf_comet_inner_res <- inner_join(msgf, comet, by="x")
msgf_comet_inner_res$error <- abs(msgf_comet_inner_res$y.x - msgf_comet_inner_res$y.y)
msgf_comet_inner_res <- msgf_comet_inner_res %>%
  filter(error <= 5)
msgf_comet_inner_res$y <- (msgf_comet_inner_res$y.x + msgf_comet_inner_res$y.y)/2
msgf_comet_inner_res <- msgf_comet_inner_res %>%
  select(x, y)
#### Three softwares overlap
all_overlap <- inner_join(msgf_xtandem_inner_res, comet_xtandem_inner_res, by="x")
all_overlap$error <- abs(all_overlap$y.x - all_overlap$y.y)
all_overlap_5 <- all_overlap %>%
  filter(error <= 5)
all_overlap_5$y <- (all_overlap_5$y.x + all_overlap_5$y.y)/2
all_overlap_5 <- all_overlap_5 %>%
  select(x, y)
msgf_xtandem_anti <- anti_join(msgf_xtandem_inner_res, all_overlap, by="x")
comet_xtandem_anti <- anti_join(comet_xtandem_inner_res, all_overlap, by="x")
msgf_comet_anti <- anti_join(msgf_comet_inner_res, all_overlap, by="x")
#### Final result
final_result <- bind_rows(all_overlap_5, msgf_xtandem_anti, comet_xtandem_anti, msgf_comet_anti)
write.table(final_result, out_file, row.names=FALSE, quote=FALSE, sep="\t")
