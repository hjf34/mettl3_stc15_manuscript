library(tidyverse)
library(ggplot2)

spike_in_table = read.table("processing/glori/resources/spike_in/spike_in_table.txt", header = T)
row.names(spike_in_table) = paste(spike_in_table$SI, spike_in_table$pos, sep = "_")

################################################################################

sample_info = as_tibble(read.csv("processing/glori/sample_info/glori_sample_info.csv"))

################################################################################

setwd("mpileup/spikein")
lf = list.files(pattern = "all_mpileup")
all_mpileup_list = lapply(lf, read.table)

acgt_counter = function(ax_mpileup){
  base_count = lapply(strsplit(ax_mpileup$V6, ","), as.numeric)
  base_names = strsplit(paste(ax_mpileup$V3, ax_mpileup$V4, sep = ","), ",")
  acgt_count = t(mapply(function(X,Y) {names(X) = Y; X[c("A","C","G","T")]}, X=base_count, Y=base_names))
  acgt_count[is.na(acgt_count)] = 0
  colnames(acgt_count) = c("A","C","G","T")
  row.names(acgt_count) = paste(ax_mpileup$V1, ax_mpileup$V2, sep = "_")
  return(acgt_count)
}

si_acgt_list = lapply(all_mpileup_list, function(n) cbind(spike_in_table, acgt_counter(n)))
names(si_acgt_list) = sample_names
si_acgt_df = do.call(rbind, si_acgt_list)
si_acgt_df$code = factor(gsub("[.].*", "", row.names(si_acgt_df)))
m6a_proportion_spike_in_table = as_tibble(si_acgt_df) %>% mutate(ACGT_count = A + C + G + T, AG_count = A + G, Arate = A/AG_count) %>% 
  left_join(sample_info, by = c("code")) %>% 
  mutate(treatx = factor(ifelse(grepl("DMSO", treat), "DMSO", "1.0 \U00B5M\nSTM3675"), levels = c("DMSO", "1.0 \U00B5M\nSTM3675"))) %>%
  mutate(repx = factor(rep)) %>%
  arrange(cell, treatx)

save(m6a_proportion_spike_in_table, file = "processing/glori/m6a_proportion_tables/m6a_proportion_spike_in_table.Rda")

################################################################################