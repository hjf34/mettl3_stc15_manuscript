library(DESeq2)
library(tidyverse)
library(tximport)
library(GenomicRanges)
library(rtracklayer)

################################################################################

gtf = rtracklayer::import("Homo_sapiens.GRCh38.108.gtf")
i2n = as_tibble(mcols(gtf[gtf$type == "transcript"])[,c("gene_id","gene_name")]) %>% distinct()

multi_cell_sample_info = as_tibble(read.csv("processing/rnaseq/multi_cell_line_stm3675/sample_info/multicell_sample_info.csv"))
cell_lines = multi_cell_sample_info %>% pull(cell) %>% unique()

################################################################################
#####Load CAOV3 runon annotation data

caov3_joined_runon_table = as_tibble(read.table("multi_cell/feature_counts/runon/featurecounts_mane_runon_joined_chunk_above_caov3_nc10", skip=1, header = T))
names1 = gsub(".umidedup.bam", "", grep(".umidedup.bam", names(caov3_joined_runon_table), value = T))
names(caov3_joined_runon_table) = c("gene_id", "chr", "start", "end", "strand", "width", names1)

caov3_joined_runon = caov3_joined_runon_table %>% mutate(gene_id = paste0(gene_id, "_ro")) %>% dplyr::select(c(1,7:12))

################################################################################
#####Load all cell runon annotation data

allcell_joined_runon_table = as_tibble(read.table("multi_cell/feature_counts/runon/featurecounts_mane_runon_joined_chunk_above_allcells_nc10", skip=1, header = T))
names1 = gsub(".umidedup.bam", "", grep(".umidedup.bam", names(allcell_joined_runon_table), value = T))
names(allcell_joined_runon_table) = c("gene_id", "chr", "start", "end", "strand", "width", names1)

allcell_joined_runon = allcell_joined_runon_table %>% mutate(gene_id = paste0(gene_id, "_ro")) %>% dplyr::select(c(1,7:42))

################################################################################

allcell_unnormalized_count_table = load("allcell_unnormalized_count_table.Rda")
caov3_unnormalized_count_table = allcell_unnormalized_count_table %>% dplyr::select(c(1, grep("CD|CT", names(.))))

caov3_ro_and_gene = bind_rows(caov3_joined_runon, caov3_unnormalized_count_table)
allcell_ro_and_gene = bind_rows(allcell_joined_runon, allcell_unnormalized_count_table)

################################################################################
################################################################################
### CAOV3 runon differential expression analysis

cts = caov3_ro_and_gene %>% dplyr::select(2:length(.)) %>% data.frame()
row.names(cts) =  caov3_ro_and_gene$gene_id

caov3_sample_info = multi_cell_sample_info %>% filter(cell == "CAOV3")
caov3_sample_df = caov3_sample_info %>% mutate(condition = factor(treat, levels = c("DMSO", "STM3675"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")

dds = DESeqDataSetFromMatrix(countData = cts, colData = caov3_sample_df, design = ~condition)
dds = DESeq(dds)
res = as_tibble(results(dds), rownames = "gene_id") %>%
  mutate(type = ifelse(grepl("_ro", gene_id), "ro", "gene")) %>%
  mutate(gene_id = gsub("_ro", "", gene_id)) %>%
  left_join(i2n) %>% relocate(c(gene_name, type), .after = "gene_id")

caov3_res_ro = res %>% filter(type == "ro") %>% 
  dplyr::select(c(gene_id, gene_name, baseMean, log2FoldChange, padj)) %>%
  dplyr::rename(c(bm = baseMean, l2fc = log2FoldChange)) %>%
  left_join(res %>% filter(type == "gene") %>% 
              dplyr::select(c(gene_id, baseMean, log2FoldChange, padj)) %>%
              dplyr::rename(c(gene_bm = baseMean, gene_l2fc = log2FoldChange, gene_padj = padj))) %>%
  left_join(caov3_joined_runon_table %>% dplyr::select(c(gene_id, chr, strand, start, end, width))) %>%
  relocate(c(chr, strand, start, end, width), .after = gene_name) %>% 
  mutate(bmbyw = bm/width)

caov3_runon_results_table = caov3_res_ro

save(caov3_runon_results_table, file = "processing/rnaseq/multi_cell_line_stm3675/caov3_runon_results_table.Rda")

################################################################################
################################################################################
### All cell runon differential expression analysis

cts = allcell_ro_and_gene %>% dplyr::select(2:length(.)) %>% data.frame()
row.names(cts) = allcell_ro_and_gene$gene_id

################################################################################
################################################################################

res_list = list()
for(a1 in 1:length(cell_lines)){
  ##############################################################################
  cell_line = cell_lines[a1]
  cell_sample_info = multi_cell_sample_info %>% filter(cell == cell_line)
  cell_sample_df = cell_sample_info %>% mutate(condition = factor(treat, levels = c("DMSO", "STM3675"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  cts1 = cts[,row.names(cell_sample_df)]
  dds = DESeqDataSetFromMatrix(countData = cts1, colData = cell_sample_df, design = ~condition)
  dds = DESeq(dds)
  ##############################################################################
  res = as_tibble(results(dds), rownames = "gene_id") %>%
    mutate(type = ifelse(grepl("_ro", gene_id), "ro", "gene")) %>%
    mutate(gene_id = gsub("_ro", "", gene_id)) %>%
    left_join(i2n) %>% relocate(c(gene_name, type), .after = "gene_id")
  ##############################################################################
  res_ro_cell = res %>% filter(type == "ro") %>% 
    dplyr::select(c(gene_id, gene_name, baseMean, log2FoldChange, padj)) %>%
    dplyr::rename(c(bm = baseMean, l2fc = log2FoldChange)) %>%
    left_join(res %>% filter(type == "gene") %>% 
                dplyr::select(c(gene_id, baseMean, log2FoldChange, padj)) %>%
                dplyr::rename(c(gene_bm = baseMean, gene_l2fc = log2FoldChange, gene_padj = padj))) %>%
    left_join(allcell_joined_runon_table %>% dplyr::select(c(gene_id, chr, strand, start, end, width))) %>%
    relocate(c(chr, strand, start, end, width), .after = gene_name) %>% 
    mutate(bmbyw = bm/width, cell = cell_line)
  ##############################################################################
  res_list[[a1]] = res_ro_cell
}

allcell_res_ro = do.call(rbind, res_list)

allcell_runon_results_table = allcell_res_ro

save(allcell_runon_results_table, file = "processing/rnaseq/multi_cell_line_stm3675/allcell_runon_results_table.Rda")

################################################################################
################################################################################
#### dose response runon

wd = "dose_response"

dose_response_sample_info = as_tibble(read.csv("processing/rnaseq/caov3_stc15_doseresponse/sample_info/doseresponse_sample_info.csv")) %>% 
  mutate(dose_factor = factor(rep(c("0","3", "10", "30", "100","300", "1000"), 3), levels = c("0","3","10","30","100","300","1000"))) 

load("doseresponse_unnormalized_count_table.Rda")

runon_joined_table_corall_caov3 = as_tibble(read.table("dose_response/feature_counts/runon/featurecounts_mane_runon_joined_chunk1000_above_caov3_nc10", skip=1, header = T))
names1 = gsub(".*bam.|.Aligned.*", "", names(runon_joined_table_corall_caov3)[7:length(runon_joined_table_corall_caov3)])
names(runon_joined_table_corall_caov3) =  c("gene_id", "chr", "start", "end", "strand", "width", names1)

runon_joined_gene_table_corall_caov3 = runon_joined_table_corall_caov3 %>% dplyr::select(c(1, 7:length(.))) %>% 
  bind_rows(doseresponse_unnormalized_count_table %>% dplyr::select(-gene_name) %>% mutate(gene_id = paste(gene_id, "gene", sep = "_")))
cts_corall_caov3 = runon_joined_gene_table_corall_caov3 %>% dplyr::select(-gene_id) %>% data.frame()
row.names(cts_corall_caov3) = runon_joined_gene_table_corall_caov3$gene_id

ro_dds_corall_caov3 = DESeqDataSetFromMatrix(countData = cts_corall_caov3, 
                                             colData = dose_response_sample_info %>% dplyr::select(code, dose_factor, rep) %>% column_to_rownames("code"), 
                                             design = ~dose_factor)

ro_dds_corall_caov3 = DESeq(ro_dds_corall_caov3)

ro_res_3_cc <- results(ro_dds_corall_caov3, contrast = c("dose_factor", "3", "0"))
ro_res_10_cc <- results(ro_dds_corall_caov3, contrast = c("dose_factor", "10", "0"))
ro_res_30_cc <- results(ro_dds_corall_caov3, contrast = c("dose_factor", "30", "0"))
ro_res_100_cc <- results(ro_dds_corall_caov3, contrast = c("dose_factor", "100", "0"))
ro_res_300_cc <- results(ro_dds_corall_caov3, contrast = c("dose_factor", "300", "0"))
ro_res_1000_cc <- results(ro_dds_corall_caov3, contrast = c("dose_factor", "1000", "0"))

################################################################################

dose_v_0 = c(3, 10, 30, 100, 300, 1000)
ro_res_pairwise_list_cc = list(ro_res_3_cc, ro_res_10_cc, ro_res_30_cc, ro_res_100_cc, ro_res_300_cc, ro_res_1000_cc)
ro_res_pairwise_list1_cc = list()
for(a1 in 1:length(ro_res_pairwise_list_cc)){
  ##############################################################################
  gene_res_pairwise = ro_res_pairwise_list_cc[[a1]] %>% as_tibble(rownames = "gene_id") %>% filter(grepl("_gene", gene_id)) %>% dplyr::select(gene_id, log2FoldChange, padj) %>%
    mutate(dose_v_0 = dose_v_0[a1]) %>% dplyr::rename(c(gene_l2fc = log2FoldChange, gene_padj = padj)) %>% mutate(gene_id = gsub("_gene", "", gene_id))
  ###############################################################################
  ro_res_pairwise_list1_cc[[a1]] = ro_res_pairwise_list_cc[[a1]] %>% as_tibble(rownames = "gene_id") %>% filter(!grepl("_gene", gene_id)) %>%
    mutate(dose_v_0 = dose_v_0[a1]) %>% left_join(runon_joined_table_corall_caov3 %>% dplyr::select(gene_id, chr, start, end, strand)) %>% 
    left_join(i2n) %>%
    relocate(c(gene_name, chr, strand, start, end), .after = gene_id) %>% 
    left_join(gene_res_pairwise)
}

dose_response_runon_results_table = do.call(rbind, ro_res_pairwise_list1_cc)
dose_response_runon_results_table$dose_v_0 = factor(dose_response_runon_results_table$dose_v_0, levels = dose_v_0)
save(dose_response_runon_results_table, file = "processing/rnaseq/caov3_stc15_doseresponse/dose_response_runon_results_table.Rda")
