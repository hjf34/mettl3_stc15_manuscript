
library(DESeq2)
library(tximport)

gtf = rtracklayer::import("/home/harry.fischl/Documents/annotations/GRCh38/Ensembl/Homo_sapiens.GRCh38.108.gtf")
i2n = as_tibble(mcols(gtf[gtf$type == "transcript"])[,c("gene_id","gene_name")]) %>% distinct()

################################################################################

wd = "multi_cell"
salmon_files = list.files(path=paste0(wd, "/salmon_counts/"), full.names = T, pattern="quant.sf$", recursive = T)
sample_names = as.vector(do.call(rbind, strsplit(salmon_files, "salmon_counts//|/quant.sf"))[,2])
names(salmon_files) = sample_names

multi_cell_sample_info = as_tibble(read.csv("processing/rnaseq/multi_cell_line_stm3675/sample_info/multicell_sample_info.csv"))
tx2gene = read.table("processing/rnaseq/genome_setup/GRCh38/transcripts_to_genes.txt")

txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)

################################################################################

allcell_unnormalized_count_table = as_tibble(round(txi[[2]]), rownames = "gene_id") %>% left_join(i2n %>% dplyr::select(gene_id, gene_name), by = join_by(gene_id)) %>% relocate(gene_name, .after = "gene_id")
#save(allcell_unnormalized_count_table, file="allcell_unnormalized_count_table.Rda")

################################################################################

cell_lines = multi_cell_sample_info$cell %>% unique %>% sort

################################################################################

res_shrink_list = list()
for(a1 in 1:length(cell_lines)){
  cell_line = cell_lines[a1]
  cell_sample_info = multi_cell_sample_info %>% filter(cell == cell_line)
  cell_sample_df = cell_sample_info %>% mutate(condition = factor(treat, levels = c("DMSO", "STM3675"))) %>% dplyr::select(code, condition) %>% column_to_rownames("code")
  ##############################################################################
  cell_txi = lapply(txi[1:3], function(n) n[,cell_sample_info %>% pull(code)])
  cell_txi[[4]] = txi[[4]]
  names(cell_txi) = names(txi)
  ##############################################################################
  cell_dds <- DESeqDataSetFromTximport(cell_txi, colData = cell_sample_df, design = ~condition)
  cell_dds <- DESeq(cell_dds)
  ##############################################################################
  res_shrink_list[[a1]] = as_tibble(lfcShrink(cell_dds, coef = "condition_STM3675_vs_DMSO", type="apeglm"), rownames = "gene_id") %>% 
    left_join(i2n %>% dplyr::select(gene_id, gene_name), by = join_by(gene_id)) %>% relocate(gene_name, .after = "gene_id") %>%
    mutate(cell = cell_line)
}

multi_cell_gene_expression_results_shrink_table = do.call(rbind, res_shrink_list) %>% mutate(cell = factor(cell, levels = cell_lines))
save(multi_cell_gene_expression_results_shrink_table, file = "processing/rnaseq/multi_cell_line_stm3675/multi_cell_gene_expression_results_shrink_table.Rda")

################################################################################
################################################################################

wd = "dose_response"
salmon_files = list.files(path=paste0(wd, "/salmon_counts/"), full.names = T, pattern="quant.sf$", recursive = T)
sample_names = as.vector(do.call(rbind, strsplit(salmon_files, "salmon_counts//|/quant.sf"))[,2])
names(salmon_files) = sample_names

dose_response_sample_info = as_tibble(read.csv("processing/rnaseq/caov3_stc15_doseresponse/sample_info/doseresponse_sample_info.csv")) %>% 
  mutate(dose_factor = factor(rep(c("0","3", "10", "30", "100","300", "1000"), 3), levels = c("0","3","10","30","100","300","1000"))) 
tx2gene = read.table("processing/rnaseq/genome_setup/GRCh38/transcripts_to_genes.txt")

txi <- tximport(salmon_files, type = "salmon", tx2gene = tx2gene)

doses = c(dose_response_sample_info$treat %>% unique)[-1]

dds_all_samples <- DESeqDataSetFromTximport(txi, 
                                            colData = dose_response_sample_info %>% dplyr::select(code, dose_factor, rep) %>% column_to_rownames("code"), 
                                            design = ~dose_factor)
dds_all_samples <- DESeq(dds_all_samples)

doseresponse_unnormalized_count_table = as_tibble(DESeq2::counts(dds_all_samples, normalized = F), rownames = "gene_id") %>% 
  left_join(i2n %>% dplyr::select(c(gene_id, gene_name))) %>% relocate(gene_name, .after = "gene_id")

#save(doseresponse_unnormalized_count_table, file = "doseresponse_unnormalized_count_table.Rda")

################################################################################

res_3_shrink <- lfcShrink(dds_all_samples, coef = c("dose_factor_3_vs_0"))
res_10_shrink <- lfcShrink(dds_all_samples, coef = c("dose_factor_10_vs_0"))
res_30_shrink <- lfcShrink(dds_all_samples, coef = c("dose_factor_30_vs_0"))
res_100_shrink <- lfcShrink(dds_all_samples, coef = c("dose_factor_100_vs_0"))
res_300_shrink <- lfcShrink(dds_all_samples, coef = c("dose_factor_300_vs_0"))
res_1000_shrink <- lfcShrink(dds_all_samples, coef = c("dose_factor_1000_vs_0"))

################################################################################

dose_v_0 = c(3, 10, 30, 100, 300, 1000)
res_pairwise_shrink_list = list(res_3_shrink, res_10_shrink, res_30_shrink, res_100_shrink, res_300_shrink, res_1000_shrink)
res_pairwise_shrink_list1 = list()
for(a1 in 1:length(res_pairwise_shrink_list)){
  res_pairwise_shrink_list1[[a1]] = res_pairwise_shrink_list[[a1]] %>% as_tibble(rownames = "gene_id") %>% 
    left_join(i2n %>% dplyr::select(c(gene_id, gene_name))) %>% relocate(gene_name, .after = gene_id) %>%
    mutate(dose_v_0 = dose_v_0[a1])
}

res_pairwise_shrink = do.call(rbind, res_pairwise_shrink_list1)
res_pairwise_shrink$dose_v_0 = factor(res_pairwise_shrink$dose_v_0, levels = dose_v_0)

dose_response_gene_expression_results_shrink_table = res_pairwise_shrink
save(dose_response_gene_expression_results_shrink_table, file = "processing/rnaseq/caov3_stc15_doseresponse/dose_response_gene_expression_results_shrink_table.Rda")

################################################################################