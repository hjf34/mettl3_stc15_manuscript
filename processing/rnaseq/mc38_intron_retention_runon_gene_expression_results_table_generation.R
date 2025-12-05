################################################################################
################################################################################
##### IR/RO
################################################################################
#### DESeq2 and DEXseq

library(DEXSeq)
library(DESeq2)
library(tidyverse)
library(tximport)
library(GenomicRanges)

################################################################################
#### IR

load("mmus_ensembl_canonical_protein_coding_introns_with_updn_exons.Rda")

mc38_sample_info = as_tibble(read.csv("processing/rnaseq/mc38_invivo_model_stc15_timeseries/sample_info/mc38invivo_sample_info.csv"))
mc38_time_points = c(1,4,12,24,48,120)

featurecounts_introns_mc38_table = as_tibble(read.table("mc38_invivo/feature_counts/intron_retention/featurecounts_introns_mc38_invivo", header = T))
names1 = gsub(".*bam.|.Aligned.*", "", names(featurecounts_introns_mc38_table)[7:length(featurecounts_introns_mc38_table)])
names(featurecounts_introns_mc38_table) =  c("gene_id", "chr", "start", "end", "strand", "width", names1)

mc38_intron_fc_table = left_join(mmus_ensembl_canonical_protein_coding_introns_with_updn_exons, featurecounts_introns_mc38_table) %>% 
  dplyr::select(gene_id, feature_id, chr, start, end, strand, mc38_sample_info$code) %>% 
  mutate(intron_id = paste(gene_id, feature_id, sep = ":")) %>% relocate(intron_id, .after = feature_id)

blood_res_list = list()
tumour_res_list = list()
for(a1 in 1:length(mc38_time_points)){
  ##############################################################################
  time_point = mc38_time_points[a1]
  mc38_time_point_sample_info = mc38_sample_info %>% filter(time == time_point | treat == "vehicle")
  blood_time_point_sample_info = mc38_time_point_sample_info %>% filter(cell == "blood")
  tumour_time_point_sample_info = mc38_time_point_sample_info %>% filter(cell == "tumour")
  ##############################################################################
  blood_time_point_sample_df = blood_time_point_sample_info %>% mutate(condition = factor(treat, levels = c("vehicle", "stc15"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  tumour_time_point_sample_df = tumour_time_point_sample_info %>% mutate(condition = factor(treat, levels = c("vehicle", "stc15"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  blood_introns_fc_table = mc38_intron_fc_table %>% dplyr::select(gene_id, feature_id, intron_id, chr, start, end, strand, blood_time_point_sample_info$code)
  tumour_introns_fc_table = mc38_intron_fc_table %>% dplyr::select(gene_id, feature_id, intron_id, chr, start, end, strand, tumour_time_point_sample_info$code)
  ##############################################################################
  blood_introns_fc_table1 = blood_introns_fc_table %>% filter(gene_id %in% c(blood_introns_fc_table %>% dplyr::count(gene_id) %>% filter(n > 1) %>% pull(gene_id)))
  blood_cts1 = blood_introns_fc_table1 %>% dplyr::select(blood_time_point_sample_info$code) %>% as.matrix()
  row.names(blood_cts1) = blood_introns_fc_table1$intron_id
  feature_ranges = blood_introns_fc_table1 %>% makeGRangesFromDataFrame()
  names(feature_ranges) = blood_introns_fc_table1$intron_id
  ##############################################################################
  design <- formula(~sample + exon + condition:exon)
  dxd = DEXSeqDataSet(countData = blood_cts1, sampleData = blood_time_point_sample_df, design = design, 
                      featureID = as.vector(blood_introns_fc_table1$feature_id), 
                      groupID = blood_introns_fc_table1$gene_id,
                      featureRanges = feature_ranges)
  dxr = DEXSeq(dxd)
  dexseq_results = as_tibble(dxr) %>% dplyr::select(groupID, featureID, exonBaseMean, log2fold_stc15_vehicle, padj)
  names(dexseq_results) = c("gene_id", "feature_id", "ebm", "l2fc", "padj")
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = blood_cts1, colData = blood_time_point_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "intron_id") %>% dplyr::select(intron_id, log2FoldChange, padj)
  names(deseq_results) = c("intron_id", "l2fc_de", "padj_de")
  ##############################################################################
  blood_res = blood_introns_fc_table1 %>% dplyr::select(c(gene_id, feature_id, intron_id)) %>%
    left_join(mmus_ensembl_canonical_protein_coding_introns_with_updn_exons) %>% 
    left_join(dexseq_results) %>% left_join(deseq_results) %>%
    mutate(stc15_time_v_vehicle = time_point) %>%
    mutate(cell = "blood") %>%
    mutate(width = end - start + 1)
  ##############################################################################
  blood_res_list[[a1]] = blood_res
  ##############################################################################
  tumour_introns_fc_table1 = tumour_introns_fc_table %>% filter(gene_id %in% c(tumour_introns_fc_table %>% dplyr::count(gene_id) %>% filter(n > 1) %>% pull(gene_id)))
  tumour_cts1 = tumour_introns_fc_table1 %>% dplyr::select(tumour_time_point_sample_info$code) %>% as.matrix()
  row.names(tumour_cts1) = tumour_introns_fc_table1$intron_id
  feature_ranges = tumour_introns_fc_table1 %>% makeGRangesFromDataFrame()
  names(feature_ranges) = tumour_introns_fc_table1$intron_id
  ##############################################################################
  design <- formula(~sample + exon + condition:exon)
  dxd = DEXSeqDataSet(countData = tumour_cts1, sampleData = tumour_time_point_sample_df, design = design, 
                      featureID = as.vector(tumour_introns_fc_table1$feature_id), 
                      groupID = tumour_introns_fc_table1$gene_id,
                      featureRanges = feature_ranges)
  dxr = DEXSeq(dxd)
  dexseq_results = as_tibble(dxr) %>% dplyr::select(groupID, featureID, exonBaseMean, log2fold_stc15_vehicle, padj)
  names(dexseq_results) = c("gene_id", "feature_id", "ebm", "l2fc", "padj")
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = tumour_cts1, colData = tumour_time_point_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "intron_id") %>% dplyr::select(intron_id, log2FoldChange, padj)
  names(deseq_results) = c("intron_id", "l2fc_de", "padj_de")
  ##############################################################################
  tumour_res = tumour_introns_fc_table1 %>% dplyr::select(c(gene_id, feature_id, intron_id)) %>%
    left_join(mmus_ensembl_canonical_protein_coding_introns_with_updn_exons) %>% 
    left_join(dexseq_results) %>% left_join(deseq_results) %>%
    mutate(stc15_time_v_vehicle = time_point) %>% 
    mutate(cell = "tumour") %>%
    mutate(width = end - start + 1)
  ##############################################################################
  tumour_res_list[[a1]] = tumour_res 
}

mc38_invivo_intron_retention_results_table = bind_rows(do.call(rbind, blood_res_list), do.call(rbind, tumour_res_list)) %>%
  dplyr::select(gene_id, gene_name, feature_id, intron_id, chr, start, end, width, strand,
                ue_range, ue_type, ue_width, de_range, de_type, de_width, 
                ebm, padj, l2fc, l2fc_de, padj_de, stc15_time_v_vehicle, cell) %>%
  dplyr::rename(c(feature = feature_id, feature_id = intron_id))

#save(mc38_invivo_intron_retention_results_table, file = "processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_intron_retention_results_table.Rda")
mc38_invivo_blood_intron_retention_results_table = mc38_invivo_intron_retention_results_table %>% filter(cell == "blood")
mc38_invivo_tumour_intron_retention_results_table = mc38_invivo_intron_retention_results_table %>% filter(cell == "tumour")
save(mc38_invivo_blood_intron_retention_results_table, file = "processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_blood_intron_retention_results_table.Rda")
save(mc38_invivo_tumour_intron_retention_results_table, file = "processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_tumour_intron_retention_results_table.Rda")

################################################################################
################################################################################
#### Gene expression analysis

mouse_gtf = rtracklayer::import("Mus_musculus.GRCm39.111.gtf")
mouse_i2n = mcols(mouse_gtf[mouse_gtf$type == "transcript"])[,c("gene_id","gene_name")]
mouse_i2n = as_tibble(mouse_i2n[!duplicated(mouse_i2n$gene_id),])
tx2gene = read.table("processing/rnaseq/genome_setup/GRCm39/transcripts_to_genes.txt")

wd = "mc38_invivo"

salmon_files = list.files(path=paste0(wd, "/salmon_counts/"), full.names = T, pattern="quant.sf$", recursive = T)
sample_names = as.vector(do.call(rbind, strsplit(salmon_files, "salmon_counts//|/quant.sf"))[,2])
names(salmon_files) = sample_names

txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)
mc38_time_points

mc38_unnormalized_count_table = as_tibble(round(txi[[2]]), rownames = "gene_id") %>% left_join(mouse_i2n %>% dplyr::select(gene_id, gene_name), by = join_by(gene_id)) %>% relocate(gene_name, .after = "gene_id")
#save(mc38_unnormalized_count_table, file="mc38_unnormalized_count_table.Rda")

################################################################################
################################################################################

blood_res_shrink_list = list()
tumour_res_shrink_list = list()
for(a1 in 1:length(mc38_time_points)){
  time_point = mc38_time_points[a1]
  mc38_time_point_sample_info = mc38_sample_info %>% filter(time == time_point | treat == "vehicle")
  blood_time_point_sample_info = mc38_time_point_sample_info %>% filter(cell == "blood")
  tumour_time_point_sample_info = mc38_time_point_sample_info %>% filter(cell == "tumour")
  ##############################################################################
  blood_time_point_sample_df = blood_time_point_sample_info %>% mutate(condition = factor(treat, levels = c("vehicle", "stc15"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  tumour_time_point_sample_df = tumour_time_point_sample_info %>% mutate(condition = factor(treat, levels = c("vehicle", "stc15"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  blood_txi = lapply(txi[1:3], function(n) n[, blood_time_point_sample_info %>% pull(code)])
  blood_txi[[4]] = txi[[4]]
  names(blood_txi) = names(txi)
  ##############################################################################
  blood_dds <- DESeqDataSetFromTximport(blood_txi, colData = blood_time_point_sample_df, design = ~condition)
  blood_dds <- DESeq(blood_dds)
  ##############################################################################
  blood_res_shrink_list[[a1]] = as_tibble(lfcShrink(blood_dds, coef = "condition_stc15_vs_vehicle", type="apeglm"), rownames = "gene_id") %>% 
    left_join(mouse_i2n %>% dplyr::select(gene_id, gene_name), by = join_by(gene_id)) %>% relocate(gene_name, .after = "gene_id") %>%
    mutate(stc15_time_v_vehicle = time_point, cell = "blood")
  ##############################################################################
  tumour_txi = lapply(txi[1:3], function(n) n[, tumour_time_point_sample_info %>% pull(code)])
  tumour_txi[[4]] = txi[[4]]
  names(tumour_txi) = names(txi)
  ##############################################################################
  tumour_dds <- DESeqDataSetFromTximport(tumour_txi, colData = tumour_time_point_sample_df, design = ~condition)
  tumour_dds <- DESeq(tumour_dds)
  ##############################################################################
  tumour_res_shrink_list[[a1]] = as_tibble(lfcShrink(tumour_dds, coef = "condition_stc15_vs_vehicle", type="apeglm"), rownames = "gene_id") %>% 
    left_join(mouse_i2n %>% dplyr::select(gene_id, gene_name), by = join_by(gene_id)) %>% relocate(gene_name, .after = "gene_id") %>%
    mutate(stc15_time_v_vehicle = time_point, cell = "tumour")
}

mc38_invivo_gene_expression_results_shrink_table = bind_rows(do.call(rbind, blood_res_shrink_list), do.call(rbind, tumour_res_shrink_list))
save(mc38_invivo_gene_expression_results_shrink_table, file = "processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_gene_expression_results_shrink_table.Rda")

################################################################################
################################################################################
################################################################################
#### RO

featurecounts_runon_mc38_table = as_tibble(read.table("mc38_invivo/feature_counts/runon/featurecounts_runon_mc38_invivo", header = T))
names1 = gsub(".*bam.|.Aligned.*", "", names(featurecounts_runon_mc38_table)[7:length(featurecounts_runon_mc38_table)])
names(featurecounts_runon_mc38_table) =  c("gene_id", "chr", "start", "end", "strand", "width", names1)

### Add gene expression counts
mc38_unnormalized_gene_count_table = mc38_unnormalized_count_table %>% mutate(gene_id = paste0(gene_id, "_gene")) %>% dplyr::select(-gene_name)

ro_cts_table = bind_rows(featurecounts_runon_mc38_table %>% dplyr::select(-c(2,3,4,5,6)), mc38_unnormalized_gene_count_table)

################################################################################

blood_ro_res_list = list()
tumour_ro_res_list = list()
for(a1 in 1:length(mc38_time_points)){
  ##############################################################################
  time_point = mc38_time_points[a1]
  mc38_time_point_sample_info = mc38_sample_info %>% filter(time == time_point | treat == "vehicle")
  blood_time_point_sample_info = mc38_time_point_sample_info %>% filter(cell == "blood")
  tumour_time_point_sample_info = mc38_time_point_sample_info %>% filter(cell == "tumour")
  ##############################################################################
  blood_time_point_sample_df = blood_time_point_sample_info %>% mutate(condition = factor(treat, levels = c("vehicle", "stc15"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  tumour_time_point_sample_df = tumour_time_point_sample_info %>% mutate(condition = factor(treat, levels = c("vehicle", "stc15"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  blood_ro_cts_table = ro_cts_table %>% dplyr::select(gene_id, blood_time_point_sample_info$code)
  tumour_ro_cts_table = ro_cts_table %>% dplyr::select(gene_id, tumour_time_point_sample_info$code)
  ##############################################################################
  blood_ro_cts = blood_ro_cts_table %>% dplyr::select(-gene_id) %>% as.matrix()
  row.names(blood_ro_cts) = blood_ro_cts_table$gene_id
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = blood_ro_cts, colData = blood_time_point_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "gene_id") %>% dplyr::select(gene_id, baseMean, log2FoldChange, padj)
  names(deseq_results) = c("gene_id", "bm", "l2fc", "padj")
  gene_res = deseq_results %>% filter(grepl("_gene", gene_id)) %>% mutate(gene_id = gsub("_gene", "", gene_id)) %>% dplyr::select(-bm)
  names(gene_res) = c("gene_id", "gene_l2fc", "gene_padj")
  ro_res = deseq_results %>% filter(!grepl("_gene", gene_id))
  blood_ro_res_list[[a1]] = left_join(ro_res, gene_res) %>% left_join(mouse_i2n) %>% relocate(gene_name, .after = gene_id) %>%
    mutate(stc15_time_v_vehicle = time_point) %>%
    mutate(cell = "blood")
  ##############################################################################
  blood_ro_cts = blood_ro_cts_table %>% dplyr::select(-gene_id) %>% as.matrix()
  row.names(blood_ro_cts) = blood_ro_cts_table$gene_id
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = blood_ro_cts, colData = blood_time_point_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "gene_id") %>% dplyr::select(gene_id, baseMean, log2FoldChange, padj)
  names(deseq_results) = c("gene_id", "bm", "l2fc", "padj")
  gene_res = deseq_results %>% filter(grepl("_gene", gene_id)) %>% mutate(gene_id = gsub("_gene", "", gene_id)) %>% dplyr::select(-bm)
  names(gene_res) = c("gene_id", "gene_l2fc", "gene_padj")
  ro_res = deseq_results %>% filter(!grepl("_gene", gene_id))
  blood_ro_res_list[[a1]] = left_join(ro_res, gene_res) %>% left_join(mouse_i2n) %>% relocate(gene_name, .after = gene_id) %>%
    mutate(stc15_time_v_vehicle = time_point) %>%
    mutate(cell = "blood")
  ##############################################################################
  tumour_ro_cts = tumour_ro_cts_table %>% dplyr::select(-gene_id) %>% as.matrix()
  row.names(tumour_ro_cts) = tumour_ro_cts_table$gene_id
  ##############################################################################
  dds = DESeqDataSetFromMatrix(countData = tumour_ro_cts, colData = tumour_time_point_sample_df, design = ~condition)
  dds = DESeq(dds)
  deseq_results = as_tibble(results(dds), rownames = "gene_id") %>% dplyr::select(gene_id, baseMean, log2FoldChange, padj)
  names(deseq_results) = c("gene_id", "bm", "l2fc", "padj")
  gene_res = deseq_results %>% filter(grepl("_gene", gene_id)) %>% mutate(gene_id = gsub("_gene", "", gene_id)) %>% dplyr::select(-bm)
  names(gene_res) = c("gene_id", "gene_l2fc", "gene_padj")
  ro_res = deseq_results %>% filter(!grepl("_gene", gene_id))
  tumour_ro_res_list[[a1]] = left_join(ro_res, gene_res) %>% left_join(mouse_i2n) %>% relocate(gene_name, .after = gene_id) %>%
    mutate(stc15_time_v_vehicle = time_point) %>%
    mutate(cell = "tumour")
}

mc38_invivo_runon_results_table = bind_rows(do.call(rbind, blood_ro_res_list), do.call(rbind, tumour_ro_res_list)) %>%
  arrange(cell, stc15_time_v_vehicle, gene_id) %>% left_join(featurecounts_runon_mc38_table %>% dplyr::select(c(1:6))) %>% 
  relocate(c("chr", "start", "end", "strand", "width"), .after = "gene_name")

save(mc38_invivo_runon_results_table, file = "processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_runon_results_table.Rda")
