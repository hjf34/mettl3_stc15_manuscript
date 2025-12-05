################################################################################
#### IR table generation
################################################################################
################################################################################
#### DESeq2 and DEXseq

library(DEXSeq)
library(DESeq2)
library(tidyverse)

################################################################################
################################################################################
#### Multi cell

featurecounts_introns_multicell_table = as_tibble(read.table("multi_cell/feature_counts/intron_retention/featurecounts_introns_multicell", header = T))
names1 = gsub(".umidedup.bam", "", grep(".umidedup.bam", names(featurecounts_introns_multicell_table), value = T))
names(featurecounts_introns_multicell_table) = c("gene_id", "chr", "start", "end", "strand", "width", names1)

load("hsap_mane_select_introns_with_updn_exons.Rda")

multi_cell_sample_info = as_tibble(read.csv("processing/rnaseq/multi_cell_line_stm3675/sample_info/multicell_sample_info.csv"))
cell_lines = multi_cell_sample_info %>% pull(cell) %>% unique()
cell_code = multi_cell_sample_info %>% mutate(cell_code = gsub("..$","", code)) %>% dplyr::select(cell, cell_code) %>% distinct()

introns_fc_table = left_join(hsap_mane_select_introns_with_updn_exons, featurecounts_introns_multicell_table) %>% 
  dplyr::select(gene_id, feature_id, chr, start, end, strand, multi_cell_sample_info$code) %>% 
  mutate(intron_id = paste(gene_id, feature_id, sep = ":")) %>% relocate(intron_id, .after = feature_id)

################################################################################

cts = introns_fc_table %>% dplyr::select(8:length(.)) %>% data.frame()
row.names(cts) = introns_fc_table$intron_id
sample_df = data.frame(row.names = colnames(cts), condition = factor(gsub("[[:digit:]]", "", colnames(cts))))

dds = DESeqDataSetFromMatrix(countData = cts, colData = sample_df, design = ~condition)
dds = estimateSizeFactors(dds)
intron_norm_counts = as_tibble(counts(dds, normalized = T), rownames = "intron_id")

################################################################################

inc = intron_norm_counts %>% pivot_longer(cols = c(2:length(.))) 
inc1 = inc %>% filter(grepl("D1", name)) %>% dplyr::rename(c("D1" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell)) %>% 
  left_join(inc %>% filter(grepl("D2", name)) %>% dplyr::rename(c("D2" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(inc %>% filter(grepl("D3", name)) %>% dplyr::rename(c("D3" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(inc %>% filter(grepl("T1", name)) %>% dplyr::rename(c("T1" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(inc %>% filter(grepl("T2", name)) %>% dplyr::rename(c("T2" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(inc %>% filter(grepl("T3", name)) %>% dplyr::rename(c("T3" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  dplyr::rename(cell_code = cell) %>%
  left_join(cell_code) %>% dplyr::select(-cell_code) %>%
  mutate(Tmax = pmax(T1, T2, T3)) %>% dplyr::select(intron_id, cell, Tmax)

################################################################################

res_list = list()
for(a1 in 1:length(cell_lines)){
  ##############################################################################
  cell_line = cell_lines[a1]
  cell_sample_info = multi_cell_sample_info %>% filter(cell == cell_line)
  cell_sample_df = cell_sample_info %>% mutate(condition = factor(treat, levels = c("DMSO", "STM3675"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  cell_introns_fc_table = introns_fc_table %>% dplyr::select(gene_id, feature_id, intron_id, chr, start, end, strand, cell_sample_info$code)
  cts = cell_introns_fc_table %>% dplyr::select(cell_sample_info$code) %>% as.matrix()
  row.names(cts) = cell_introns_fc_table$intron_id
  cell_introns_fc_table1 = cell_introns_fc_table %>% filter(gene_id %in% c(cell_introns_fc_table %>% dplyr::count(gene_id) %>% filter(n > 1) %>% pull(gene_id)))
  cts1 = cell_introns_fc_table1 %>% dplyr::select(cell_sample_info$code) %>% as.matrix()
  row.names(cts1) = cell_introns_fc_table1$intron_id
  feature_ranges = cell_introns_fc_table1 %>% makeGRangesFromDataFrame()
  names(feature_ranges) = cell_introns_fc_table1$intron_id
  ##############################################################################
  design <- formula(~sample + exon + condition:exon)
  dxd = DEXSeqDataSet(countData = cts1, sampleData = cell_sample_df, design =  design, 
                      featureID = as.vector(cell_introns_fc_table1$feature_id), 
                      groupID = cell_introns_fc_table1$gene_id,
                      featureRanges = feature_ranges)
  dxr = DEXSeq(dxd)
  dexseq_results = as_tibble(dxr) %>% dplyr::select(groupID, featureID, exonBaseMean, log2fold_STM3675_DMSO, padj)
  names(dexseq_results) = c("gene_id", "feature_id", "ebm", "l2fc", "padj")
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = cts, colData = cell_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "intron_id") %>% dplyr::select(intron_id, log2FoldChange, padj)
  names(deseq_results) = c("intron_id", "l2fc_de", "padj_de")
  ##############################################################################
  res = cell_introns_fc_table1 %>% dplyr::select(c(gene_id, feature_id, intron_id)) %>%
    left_join(hsap_mane_select_introns_with_updn_exons) %>% 
    left_join(dexseq_results) %>% left_join(deseq_results) %>%
    mutate(cell = cell_line) %>% left_join(inc1) %>%
    mutate(width = end - start + 1) %>%
    mutate(tbyw = Tmax/width)
  ##############################################################################
  res_list[[a1]] = res
}

allcell_intron_retention_results_table =  do.call(rbind, res_list) %>% 
  separate(ue_range, into = c("ue_chr", "ue_se", "ue_strand"), sep = ":") %>%
  separate(ue_se, into = c("ue_start", "ue_end"), sep = "-") %>% 
  separate(de_range, into = c("de_chr", "de_se", "de_strand"), sep = ":") %>%
  separate(de_se, into = c("de_start", "de_end"), sep = "-") %>%
  mutate(genomic_data = paste0("chr", chr, ":", start, "-", end, ":", strand)) %>%
  dplyr::select(cell, gene_id, gene_name, feature_id, intron_id, chr, start, end, width, strand,
                genomic_data, ue_start, ue_end, ue_width, ue_type, de_start, de_end, de_width, de_type, 
                ebm, padj, l2fc, l2fc_de, padj_de, Tmax, tbyw) %>%
  dplyr::rename(c(feature = feature_id, feature_id = intron_id))

#save(allcell_intron_retention_results_table, file = "processing/rnaseq/multi_cell_line_stm3675/allcell_intron_retention_results_table.Rda")

################################################################################
################################################################################
################################################################################
#### Dose response

featurecounts_introns_doseresponse_table = as_tibble(read.table("dose_response/feature_counts/intron_retention/featurecounts_introns_doseresponse", header = T))
names1 = gsub(".*bam.|.Aligned.*", "", names(featurecounts_introns_doseresponse_table)[7:length(featurecounts_introns_doseresponse_table)])
names(featurecounts_introns_doseresponse_table) =  c("gene_id", "chr", "start", "end", "strand", "width", names1)

dose_response_sample_info = as_tibble(read.csv("processing/rnaseq/caov3_stc15_doseresponse/sample_info/doseresponse_sample_info.csv")) %>% 
  mutate(dose_factor = factor(rep(c("0","3", "10", "30", "100","300", "1000"), 3), levels = c("0","3","10","30","100","300","1000")))

dose_v_0 = c(3, 10, 30, 100, 300, 1000)

introns_dr_fc_table = left_join(hsap_mane_select_introns_with_updn_exons, featurecounts_introns_doseresponse_table) %>% 
  dplyr::select(gene_id, feature_id, chr, start, end, strand, dose_response_sample_info$code) %>% 
  mutate(intron_id = paste(gene_id, feature_id, sep = ":")) %>% relocate(intron_id, .after = feature_id)

res_list = list()
for(a1 in 1:length(dose_v_0)){
  ##############################################################################
  dose = dose_v_0[a1]
  dose_sample_info = dose_response_sample_info %>% filter(treat %in% c(0, dose)) %>% arrange(treat)
  dose_sample_df = dose_sample_info %>% mutate(condition = factor(treat, levels = c(0, dose))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  dose_introns_fc_table = introns_dr_fc_table %>% dplyr::select(gene_id, feature_id, intron_id, cell_sample_info$code)
  cts =  dose_introns_fc_table %>% dplyr::select(dose_sample_info$code) %>% as.matrix()
  row.names(cts) = dose_introns_fc_table$intron_id
  dose_introns_fc_table1 = dose_introns_fc_table %>% filter(gene_id %in% c(dose_introns_fc_table %>% dplyr::count(gene_id) %>% filter(n > 1) %>% pull(gene_id)))
  cts1 = dose_introns_fc_table1 %>% dplyr::select(dose_sample_info$code) %>% as.matrix()
  row.names(cts1) = dose_introns_fc_table1$intron_id
  feature_ranges = dose_introns_fc_table1 %>% makeGRangesFromDataFrame()
  names(feature_ranges) = dose_introns_fc_table1$intron_id
  ##############################################################################
  design <- formula(~sample + exon + condition:exon)
  dxd = DEXSeqDataSet(countData = cts1, sampleData = dose_sample_df, design = design, 
                      featureID = as.vector(dose_introns_fc_table1$feature_id), 
                      groupID = dose_introns_fc_table1$gene_id,
                      featureRanges = feature_ranges)
  dxr = DEXSeq(dxd)
  dexseq_results = as_tibble(dxr) %>% dplyr::select(c(1,2,3,10,7))
  names(dexseq_results) = c("gene_id", "feature_id", "ebm", "l2fc", "padj")
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = cts1, colData = dose_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "intron_id") %>% dplyr::select(intron_id, log2FoldChange, padj)
  names(deseq_results) = c("intron_id", "l2fc_de", "padj_de")
  ##############################################################################
  res = cell_introns_fc_table1 %>% dplyr::select(c(gene_id, feature_id, intron_id)) %>%
    left_join(hsap_mane_select_introns_with_updn_exons) %>% 
    left_join(dexseq_results) %>% left_join(deseq_results) %>%
    mutate(dose = dose) %>%
    mutate(width = end - start + 1)
  ##############################################################################
  res_list[[a1]] = res
}

doseresponse_intron_retention_results_table = do.call(rbind, res_list) %>% 
  separate(ue_range, into = c("ue_chr", "ue_se", "ue_strand"), sep = ":") %>%
  separate(ue_se, into = c("ue_start", "ue_end"), sep = "-") %>% 
  separate(de_range, into = c("de_chr", "de_se", "de_strand"), sep = ":") %>%
  separate(de_se, into = c("de_start", "de_end"), sep = "-") %>%
  mutate(genomic_data = paste0("chr", chr, ":", start, "-", end, ":", strand)) %>%
  dplyr::select(dose, gene_id, gene_name, feature_id, intron_id, chr, start, end, width, strand,
                genomic_data, ue_start, ue_end, ue_width, ue_type, de_start, de_end, de_width, de_type, 
                ebm, padj, l2fc, l2fc_de, padj_de) %>%
  dplyr::rename(c(feature = feature_id, feature_id = intron_id))

#save(dose_response_intron_retention_results_table, file = "processing/rnaseq/caov3_stc15_doseresponse/dose_response_intron_retention_results_table.Rda")
