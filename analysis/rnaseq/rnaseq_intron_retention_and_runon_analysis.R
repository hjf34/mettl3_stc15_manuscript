library(ggpubr)

load("processing/glori/m6a_proportion_tables/m6a_proportion_table_filter_ag20_m6a25.Rda")

load("processing/rnaseq/multi_cell_line_stm3675/allcell_intron_retention_results_table.Rda")
load("processing/rnaseq/multi_cell_line_stm3675/caov3_runon_results_table.Rda")
load("processing/rnaseq/multi_cell_line_stm3675/allcell_runon_results_table.Rda")

load("processing/rnaseq/caov3_stc15_doseresponse/dose_response_intron_retention_results_table.Rda")
load("processing/rnaseq/caov3_stc15_doseresponse/dose_response_runon_results_table.Rda")

load("processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_blood_intron_retention_results_table.Rda")
load("processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_tumour_intron_retention_results_table.Rda")
load("processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_runon_results_table.Rda")

mc38_invivo_intron_retention_results_table = bind_rows(mc38_invivo_blood_intron_retention_results_table, mc38_invivo_tumour_intron_retention_results_table)

################################################################################
#m6a_proportion_table_filter_ag20_m6a25 
shared_sites_ag20_m6a25 = m6a_proportion_table_filter_ag20_m6a25 %>% filter(D_m6a >= 0.25 & D_ag >= 20) %>% 
  group_by(xpos, gene_id) %>% summarise(count = n()) %>% filter(count == 6) %>% pull(xpos) %>% unique()
m6a_proportion_table_shared_sites = m6a_proportion_table_filter_ag20_m6a25 %>% filter(xpos %in% shared_sites_ag20_m6a25)

################################################################################
#### Significant intron retention events (downstream of long exons) in CAOV3 and in any cell line
caov3_intron_retention_all_uewidth_significant = allcell_intron_retention_results_table %>% filter(cell == "CAOV3", l2fc > 0.5, padj < 0.05, l2fc_de > 0, padj_de < 0.05, width > 250, tbyw > 0.1)
caov3_intron_retention_300nt_uewidth_significant = allcell_intron_retention_results_table %>% filter(cell == "CAOV3", l2fc > 0.5, padj < 0.05, l2fc_de > 0, padj_de < 0.05, width > 250, tbyw > 0.1, ue_width >= 300)
intron_retention_any_significant = allcell_intron_retention_results_table %>% filter(l2fc > 0.5, padj < 0.05, l2fc_de > 0, padj_de < 0.05, width > 250, tbyw > 0.1, ue_width >= 300) %>% pull(genomic_data) %>% unique()
allcell_intron_retention_300nt_uewidth_significant_any = allcell_intron_retention_results_table %>% filter(genomic_data %in% intron_retention_any_significant)

data_table_caov3_intron_retention_all_uewidth_significant = caov3_intron_retention_all_uewidth_significant
names(data_table_caov3_intron_retention_all_uewidth_significant)[c(12:26)] = c("upstream_exon_start", "upstream_exon_end", "upstream_exon_width", "upstream_exon_type",
                                                                               "downstream_exon_start", "downstream_exon_end", "downstream_exon_width", "downstream_exon_type",
                                                                               "exon_base_mean", "DEXseq_padj", "DEXseq_log2FC", "DESeq_l2fc", "DESeq_padj", 
                                                                               "treated_max_count", "treated_max_count_by_intron_width")
data_table_caov3_intron_retention_300nt_uewidth_significant = caov3_intron_retention_300nt_uewidth_significant
names(data_table_caov3_intron_retention_300nt_uewidth_significant) = names(data_table_caov3_intron_retention_all_uewidth_significant)
data_table_allcell_intron_retention_300nt_uewidth_significant_any = allcell_intron_retention_300nt_uewidth_significant_any
names(data_table_allcell_intron_retention_300nt_uewidth_significant_any) = names(data_table_caov3_intron_retention_all_uewidth_significant)

write.csv(data_table_caov3_intron_retention_all_uewidth_significant, file = "analysis/rnaseq/tables/data_table_significant_intron_retention_events_in_caov3_all_upstream_exon_widths.csv", quote = F, row.names = F)
write.csv(data_table_caov3_intron_retention_300nt_uewidth_significant, file = "analysis/rnaseq/tables/data_table_significant_intron_retention_events_in_caov3_300nt_upstream_exon_widths.csv", quote = F, row.names = F)
write.csv(data_table_allcell_intron_retention_300nt_uewidth_significant_any, file = "analysis/rnaseq/tables/data_table_significant_intron_retention_events_in_any_cell_300nt_upstream_exon_widths.csv", quote = F, row.names = F)

################################################################################
#### Significant runon events in CAOV3 and in any cell line
caov3_runon_significant = caov3_runon_results_table %>% filter(padj < 0.05, l2fc > 1, gene_l2fc < 1)
runon_any_significant = allcell_runon_results_table %>% filter(padj < 0.05, l2fc > 1, gene_l2fc < 1) %>% pull(gene_id) %>% unique()
allcell_runon_significant_any = allcell_runon_results_table %>% filter(gene_id %in% runon_any_significant)

data_table_caov3_runon_significant = caov3_runon_significant
names(data_table_caov3_runon_significant)[3:14] = c("runon_chr", "runon_strand", "runon_start", "runon_end", "runon_width", "runon_base_mean", "runon_log2FC", "runon_padj", "gene_base_mean", "gene_log2FC", "gene_padj", "base_mean_by_runon_width")
data_table_allcell_runon_significant_any = allcell_runon_significant_any
names(data_table_allcell_runon_significant_any)[3:15] = c("runon_chr", "runon_strand", "runon_start", "runon_end", "runon_width", "runon_base_mean", "runon_log2FC", "runon_padj", "gene_base_mean", "gene_log2FC", "gene_padj", "base_mean_by_runon_width", "cell")

write.csv(data_table_caov3_runon_significant, file = "analysis/rnaseq/tables/data_table_significant_runon_events_in_caov3.csv", quote = F, row.names = F)
write.csv(data_table_allcell_runon_significant_any, file = "analysis/rnaseq/tables/data_table_significant_runon_events_in_any_cell.csv", quote = F, row.names = F)

caov3_intron_retention_300nt_uewidth_significant = caov3_intron_retention_300nt_uewidth_significant %>%
  mutate(exon_annotation = paste0("chr", chr, ":", ue_start, "-", ue_end, ":", strand))
allcell_intron_retention_300nt_uewidth_significant_any = allcell_intron_retention_300nt_uewidth_significant_any %>%
  mutate(exon_annotation = paste0("chr", chr, ":", ue_start, "-", ue_end, ":", strand))

#### Expected number of intron retention events downstream of long exons v actual number
exon_width_stats = allcell_intron_retention_results_table %>% filter(cell == "CAOV3", width > 250, tbyw > 0.1) %>% pull(ue_width) 
round(dim(caov3_intron_retention_all_uewidth_significant)[1]*sum(exon_width_stats >= 300)/length(exon_width_stats)) # 17 = expected
dim(caov3_intron_retention_300nt_uewidth_significant)[1] # 79 = actual

################################################################################
################################################################################
### GO:BP analysis

gobp_gene_sets =  msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()

caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result = clusterProfiler::enricher(gene = caov3_intron_retention_300nt_uewidth_significant$gene_name, TERM2GENE = gobp_gene_sets)
caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result_table = caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result %>% data.frame() %>% as_tibble() %>% dplyr::select(-Description)

caov3_runon_significant_gobp_enrich_result = clusterProfiler::enricher(gene = caov3_runon_significant$gene_name, TERM2GENE = gobp_gene_sets)
caov3_runon_significant_gobp_enrich_result_table = caov3_runon_significant_gobp_enrich_result %>% data.frame() %>% as_tibble() %>% dplyr::select(-Description)

################################################################################

save(caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result, file = "analysis/rnaseq/tables/caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result.Rda")
ggsave(plot = enrichplot::dotplot(caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result),
       filename =  "analysis/rnaseq/plots/caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result_dotplot.png",
       height = 20, width = 20, units = "cm", dpi = 300)
write.csv(caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result_table, 
          file = "analysis/rnaseq/tables/caov3_intron_retention_300nt_uewidth_significant_gobp_enrich_result_table.csv", 
          quote = F, row.names = F)

save(caov3_runon_significant_gobp_enrich_result, file = "analysis/rnaseq/tables/caov3_runon_significant_gobp_enrich_result.Rda")
ggsave(plot = enrichplot::dotplot(caov3_runon_significant_gobp_enrich_result),
       filename =  "analysis/rnaseq/plots/caov3_runon_significant_gobp_enrich_result_dotplot.png",
       height = 18, width = 18, units = "cm", dpi = 300)
write.csv(caov3_runon_significant_gobp_enrich_result_table, 
          file = "analysis/rnaseq/tables/caov3_runon_significant_gobp_enrich_result_table.csv", 
          quote = F, row.names = F)

################################################################################
################################################################################

allcell_intron_retention_300nt_uewidth_significant_any_boxplot = 
  ggplot(allcell_intron_retention_300nt_uewidth_significant_any, aes(x = cell, y = l2fc, fill = cell))+
  geom_boxplot(outlier.shape = NA, notch = T, alpha = 0.8)+
  coord_cartesian(ylim = c(0,2))+
  geom_jitter(alpha = 0.2, width = 0.1, size= 0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Cell", y = "Log2 FC (1.0 µM STM3675/DMSO) in intron retention\n(significant events in any cell line)", fill = "Cell")
  
allcell_runon_significant_any_boxplot = 
  ggplot(allcell_runon_significant_any, aes(x = cell, y = l2fc, fill = cell))+
  geom_boxplot(outlier.shape = NA, notch = T, alpha = 0.8)+
  coord_cartesian(ylim = c(0,3))+
  geom_jitter(alpha = 0.2, width = 0.1, size= 0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Cell", y = "Log2 FC (1.0 µM STM3675/DMSO) in run-on\n(significant events in any cell line)", fill = "Cell") 

################################################################################

ggsave(plot = allcell_intron_retention_300nt_uewidth_significant_any_boxplot,
       filename =  "analysis/rnaseq/plots/allcell_intron_retention_300nt_uewidth_significant_any_boxplot.png",
       height = 12.5, width = 10, units = "cm", dpi = 300)
ggsave(plot = allcell_runon_significant_any_boxplot,
       filename =  "analysis/rnaseq/plots/allcell_runon_significant_any_boxplot.png",
       height = 12.5, width = 10, units = "cm", dpi = 300)

ggsave(plot = allcell_intron_retention_300nt_uewidth_significant_any_boxplot + stat_compare_means(method = "kruskal.test", label.y = 2, label.x = 2, size = 3), # p < 2.2e-16
       filename =  "analysis/rnaseq/plots/allcell_intron_retention_300nt_uewidth_significant_any_boxplot_with_stats.png",
       height = 12.5, width = 10, units = "cm", dpi = 300)
ggsave(plot = allcell_runon_significant_any_boxplot + stat_compare_means(method = "kruskal.test", label.y = 3, label.x = 2, size = 3), # p < 2.2e-16
       filename =  "analysis/rnaseq/plots/allcell_runon_significant_any_boxplot_with_stats.png",
       height = 12.5, width = 10, units = "cm", dpi = 300)

################################################################################
################################################################################
#### m6A analysis

m6a_proportion_filter = 0.25
ag_count_filter = 20

m6a_proportion_table_filter_ag20_m6a25_caov3 = m6a_proportion_table_filter_ag20_m6a25 %>% filter(cell %in% c("CAOV3")) %>% 
  filter(D_m6a >= m6a_proportion_filter & D_ag >= ag_count_filter) %>% 
  dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>% distinct()

################################################################################
################################################################################

exon_width1 = 300
m6a_proportion_table_filter_ag20_m6a25_caov3 = m6a_proportion_table_filter_ag20_m6a25 %>% 
  filter(cell %in% c("CAOV3")) %>% 
  filter(!is.na(gene_id)) %>% 
  filter(D_m6a >= m6a_proportion_filter & D_ag >= ag_count_filter)

m6a_caov3_irro_annotation = m6a_proportion_table_filter_ag20_m6a25_caov3 %>%
  mutate(exon_type1 = paste0(ifelse(exon_type %in% c("first", "internal"), 
                                    "no IR upstream first/internal exon", 
                                    ifelse(exon_type %in% "solo", "no RO upstream single exon", "no RO upstream terminal exon")), 
                             ifelse(exon_width >= exon_width1, sprintf(" \u2265 %s nt", exon_width1), sprintf(" < %s nt", exon_width1)))) %>%
  mutate(exon_type2 = ifelse(exon_annotation %in% caov3_intron_retention_300nt_uewidth_significant$exon_annotation, 
                             sprintf("IR upstream first/internal exon \u2265 %s nt", exon_width1), 
                             ifelse(gene_id %in% caov3_runon_significant$gene_id,
                                    ifelse(exon_type %in% c("last"), 
                                           sprintf("RO upstream terminal exon \u2265 %s nt", exon_width1), 
                                           ifelse(exon_type %in% c("solo"), 
                                                  sprintf("RO upstream single exon \u2265 %s nt", exon_width1), 
                                                  as.vector(exon_type1))),
                                    as.vector(exon_type1)))) %>%
  mutate(exon_type3 = factor(exon_type2, 
                             levels = c("no IR upstream first/internal exon < 300 nt", 
                                        "no IR upstream first/internal exon ≥ 300 nt", "IR upstream first/internal exon ≥ 300 nt", 
                                        "no RO upstream terminal exon < 300 nt", "no RO upstream terminal exon ≥ 300 nt", "RO upstream terminal exon ≥ 300 nt",
                                        "no RO upstream single exon ≥ 300 nt", "RO upstream single exon ≥ 300 nt"))) %>%
  mutate(facet1 = ifelse(grepl("RO", exon_type3), "Run-on", "Intron retention"),
         col1 = ifelse(grepl("^no", exon_type3), "no", "yes"),
         col2 = factor(ifelse(grepl("^no", exon_type3), ifelse(grepl("Run-on", facet1), "no RO", "no IR"), ifelse(grepl("Run-on", facet1), "RO", "IR")), levels = c("no RO", "no IR", "RO", "IR")))

m6a_caov3_irro_annotation$exon_type4 = m6a_caov3_irro_annotation$exon_type3
levels(m6a_caov3_irro_annotation$exon_type4) = gsub("exon ", "exon\n", gsub("RO ", "RO\n", gsub("IR ", "IR\n", gsub("upstream first/internal exon |upstream ", "", levels(m6a_caov3_irro_annotation$exon_type3)))))

caov3_m6a_site_upstream_ir_boxplot = ggplot(m6a_caov3_irro_annotation %>% filter(grepl("IR", exon_type3)), aes(x = exon_type4, y = D_m6a, fill = col2))+
  geom_boxplot(notch = T, coef = 50, show.legend = T)+
  facet_grid(~facet1)+
  theme_bw()+
  coord_cartesian(ylim = c(0,1.05))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Location of m6A-modified site", y = "m6A proportion",
       fill = "IR event\ndownstream")+
  stat_compare_means(comparisons = list(c("no IR\n≥ 300 nt", "IR\n≥ 300 nt")), 
                     label.y = 1.01, method = "wilcox.test")

caov3_m6a_site_upstream_ro_boxplot = ggplot(m6a_caov3_irro_annotation %>% filter(grepl("RO", exon_type3)), aes(x = exon_type4, y = D_m6a, fill = col2))+
  geom_boxplot(notch = T, coef = 50, show.legend = T)+
  facet_grid(~facet1)+
  theme_bw()+
  coord_cartesian(ylim = c(0,1.05))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Location of m6A-modified site", y = "m6A proportion",
       fill = "RO event\ndownstream")+
  stat_compare_means(comparisons = list(c("no RO\nterminal exon\n≥ 300 nt", "RO\nterminal exon\n≥ 300 nt"),
                                        c("no RO\nsingle exon\n≥ 300 nt", "RO\nsingle exon\n≥ 300 nt")), 
                     label.y = 1.01, method = "wilcox.test")

ggsave(filename = "analysis/rnaseq/plots/caov3_m6a_site_upstream_ir_boxplot.png",
       plot = caov3_m6a_site_upstream_ir_boxplot, width = 8, height = 13, units = "cm", dpi = 300)
ggsave(filename = "analysis/rnaseq/plots/caov3_m6a_site_upstream_ro_boxplot.png",
       plot = caov3_m6a_site_upstream_ro_boxplot, width = 11, height = 13, units = "cm", dpi = 300)

################################################################################
################################################################################

shared_sites_ag20_m6a25 = m6a_proportion_table_filter_ag20_m6a25 %>% filter(D_m6a >= 0.25 & D_ag >= 20) %>% 
  group_by(xpos, gene_id) %>% summarise(count = n()) %>% filter(count == 6) %>% pull(xpos) %>% unique()
m6a_proportion_table_shared_sites = m6a_proportion_table_filter_ag20_m6a25 %>% filter(xpos %in% shared_sites_ag20_m6a25)
m6a_proportion_table_shared_sites_distinct = m6a_proportion_table_shared_sites %>% 
  dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>% distinct()

################################################################################
################################################################################

m6a_proportion_table_shared_sites_irro_annotation = m6a_proportion_table_shared_sites %>% filter(!is.na(gene_id)) %>%
  mutate(exon_type1 = paste0(ifelse(exon_type %in% c("first", "internal"), 
                                    "no IR upstream first/internal exon", 
                                    ifelse(exon_type %in% "solo", "no RO upstream single exon", "no RO upstream terminal exon")), 
                             ifelse(exon_width >= exon_width1, sprintf(" \u2265 %s nt", exon_width1), sprintf(" < %s nt", exon_width1)))) %>%
  mutate(exon_type2 = ifelse(exon_annotation %in% allcell_intron_retention_300nt_uewidth_significant_any$exon_annotation, 
                             sprintf("IR upstream first/internal exon \u2265 %s nt", exon_width1), 
                             ifelse(gene_id %in% allcell_runon_significant_any$gene_id,
                                    ifelse(exon_type %in% c("last"), 
                                           sprintf("RO upstream terminal exon \u2265 %s nt", exon_width1), 
                                           ifelse(exon_type %in% c("solo"), 
                                                  sprintf("RO upstream single exon \u2265 %s nt", exon_width1), 
                                                  as.vector(exon_type1))),
                                    as.vector(exon_type1)))) %>%
  mutate(exon_type3 = factor(exon_type2, 
                             levels = c("no IR upstream first/internal exon < 300 nt", 
                                        "no IR upstream first/internal exon ≥ 300 nt", "IR upstream first/internal exon ≥ 300 nt", 
                                        "no RO upstream terminal exon < 300 nt", "no RO upstream terminal exon ≥ 300 nt", "RO upstream terminal exon ≥ 300 nt",
                                        "no RO upstream single exon ≥ 300 nt", "RO upstream single exon ≥ 300 nt"))) %>%
  mutate(facet1 = ifelse(grepl("RO", exon_type3), "Run-on", "Intron retention"))

m6a_proportion_table_shared_sites_irro_annotation$exon_type4 = m6a_proportion_table_shared_sites_irro_annotation$exon_type3
levels(m6a_proportion_table_shared_sites_irro_annotation$exon_type4) = gsub("exon ", "exon\n", gsub("RO ", "RO\n", gsub("IR ", "IR\n", gsub("upstream first/internal exon |upstream ", "", levels(m6a_proportion_table_shared_sites_irro_annotation$exon_type3)))))

################################################################################

allcell_m6a_site_upstream_ir_boxplot = ggplot(m6a_proportion_table_shared_sites_irro_annotation %>% filter(grepl("IR", exon_type3)), 
                                              aes(x = exon_type4, y = D_m6a, fill = cell))+
  geom_boxplot(notch = T, coef = 50, show.legend = T)+
  facet_grid(~facet1)+
  theme_bw()+
  coord_cartesian(ylim = c(0,1.05))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Location of m6A-modified site", y = "m6A proportion",
       fill = "Cell")

allcell_m6a_site_upstream_ro_boxplot = ggplot(m6a_proportion_table_shared_sites_irro_annotation %>% filter(grepl("RO", exon_type3)), 
                                              aes(x = exon_type4, y = D_m6a, fill = cell))+
  geom_boxplot(notch = T, coef = 50, show.legend = T)+
  facet_grid(~facet1)+
  theme_bw()+
  coord_cartesian(ylim = c(0,1.05))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Location of m6A-modified site", y = "m6A proportion",
       fill = "Cell")

ggsave(filename = "analysis/rnaseq/plots/allcell_m6a_site_upstream_ir_boxplot.png",
       plot = allcell_m6a_site_upstream_ir_boxplot, width = 11, height = 13, units = "cm", dpi = 300)
ggsave(filename = "analysis/rnaseq/plots/allcell_m6a_site_upstream_ro_boxplot.png",
       plot = allcell_m6a_site_upstream_ro_boxplot, width = 14, height = 13, units = "cm", dpi = 300)

################################################################################

m6a_proportion_table_shared_sites_roterminalandsolo_annotation = m6a_proportion_table_shared_sites_irro_annotation %>% 
  filter(exon_type3 %in% c("RO upstream terminal exon ≥ 300 nt", "RO upstream solo exon ≥ 300 nt"))
m6a_proportion_table_shared_sites_irfi_annotation = m6a_proportion_table_shared_sites_irro_annotation %>% 
  filter(exon_type3 == "IR upstream first/internal exon ≥ 300 nt")

m6a_terminalandsolo_exon_rol2fc_correlation_plot = m6a_proportion_table_shared_sites_roterminalandsolo_annotation %>% 
  group_by(cell) %>% summarise(across(D_m6a, median)) %>% 
  left_join(allcell_runon_significant_any %>% filter(gene_id %in% m6a_proportion_table_shared_sites_roterminalandsolo_annotation$gene_id) %>%
              group_by(cell) %>% summarise(across(l2fc, median))) %>%
  ggplot(aes(x = D_m6a, y = l2fc))+
  geom_point(size = 0.1)+
  geom_smooth(method = "lm", se = F, color = "grey", lty = 2)+
  stat_cor(method="pearson")+
  geom_point(aes(col = cell), size = 3) +
  theme_bw()+
  labs(x = "Median m6A proportion at m6A-modified sites\nin exons upstream of run-on events",
       y = "Median log2FC (1.0 µM STM3675/DMSO)\nin run-on at significant events",
       col = "Cell")

m6a_upstream_exon_irl2fc_correlation_plot = m6a_proportion_table_shared_sites_irfi_annotation %>% 
  group_by(cell) %>% summarise(across(D_m6a, median)) %>% 
  left_join(allcell_intron_retention_300nt_uewidth_significant_any %>% filter(exon_annotation %in% m6a_proportion_table_shared_sites_irfi_annotation$exon_annotation) %>%
              group_by(cell) %>% summarise(across(l2fc, median))) %>%
  ggplot(aes(x = D_m6a, y = l2fc))+
  geom_point(size = 0.1)+
  geom_smooth(method = "lm", se = F, color = "grey", lty = 2)+
  stat_cor(method="pearson")+
  geom_point(aes(col = cell), size = 3) +
  theme_bw()+
  labs(x = "Median m6A proportion at m6A-modified sites\nin exons upstream of intron retention events",
       y = "Median log2FC (1.0 µM STM3675/DMSO)\nin intron retention at significant events",
       col = "Cell")

################################################################################

ggsave(filename = "analysis/rnaseq/plots/m6a_terminalandsolo_exon_rol2fc_correlation_plot.png",
       plot = m6a_terminalandsolo_exon_rol2fc_correlation_plot, width = 12, height = 12, units = "cm", dpi = 300)
ggsave(filename = "analysis/rnaseq/plots/m6a_upstream_exon_irl2fc_correlation_plot.png",
       plot = m6a_upstream_exon_irl2fc_correlation_plot, width = 12, height = 12, units = "cm", dpi = 300)

################################################################################
#### Dose response

caov3_intron_retention_300nt_uewidth_significant = allcell_intron_retention_results_table %>% 
  filter(cell == "CAOV3", l2fc > 0.5, padj < 0.05, l2fc_de > 0, padj_de < 0.05, width > 250, tbyw > 0.1, ue_width >= 300)

data_table_dose_response_intron_retention_results = dose_response_intron_retention_results_table %>% filter(genomic_data %in% caov3_intron_retention_300nt_uewidth_significant$genomic_data)
names(data_table_dose_response_intron_retention_results)[c(12:24)] = c("upstream_exon_start", "upstream_exon_end", "upstream_exon_width", "upstream_exon_type",
                                                                       "downstream_exon_start", "downstream_exon_end", "downstream_exon_width", "downstream_exon_type",
                                                                       "exon_base_mean", "DEXseq_padj", "DEXseq_log2FC", "DESeq_l2fc", "DESeq_padj")

write.csv(data_table_dose_response_intron_retention_results, file = "analysis/rnaseq/tables/data_table_significant_intron_retention_events_in_caov3_dose_response.csv", quote = F, row.names = F)

dose_response_intron_retention_boxplot = 
  ggplot(dose_response_intron_retention_results_table %>% filter(genomic_data %in% caov3_intron_retention_300nt_uewidth_significant$genomic_data), 
         aes(x = dose, y = l2fc, fill = dose))+
  geom_boxplot(outlier.shape = NA, notch = T)+
  coord_cartesian(ylim = c(-0.5,1.5))+
  geom_jitter(alpha = 0.2, width = 0.1, size= 0.5)+
  scale_fill_brewer(palette = "Greens")+
  theme_bw()+
  labs(x = "STC-15 (nM) v DMSO", y = "Log2 FC (STC-15/DMSO) in intron retention", 
       fill = "STC-15 (nM)") 

################################################################################

caov3_runon_significant = caov3_runon_results_table %>% filter(padj < 0.05, l2fc > 1, gene_l2fc < 1)

data_table_dose_response_runon_results = dose_response_runon_results_table %>% filter(gene_id %in% (caov3_runon_significant %>% pull(gene_id))) %>% dplyr::select(c(1:8,12:15))
names(data_table_dose_response_runon_results)[3:12] = c("runon_chr", "runon_strand", "runon_start", "runon_end", "runon_base_mean", "runon_log2FC", "runon_padj", "dose_v_0", "gene_log2FC", "gene_padj")

write.csv(data_table_dose_response_runon_results, file = "analysis/rnaseq/tables/data_table_significant_runon_events_in_caov3_dose_response.csv", quote = F, row.names = F)

dose_response_runon_boxplot = 
  dose_response_runon_results_table %>% filter(gene_id %in% (caov3_runon_significant %>% pull(gene_id))) %>%
  ggplot(aes(x = dose_v_0, y = log2FoldChange, fill = dose_v_0)) +
  geom_boxplot(outlier.shape = NA, notch = T) +
  scale_fill_brewer(palette = "Greens")+
  geom_jitter(alpha = 0.1, size = 0.2, col = "black", width = 0.1)+
  labs(x = "STC-15 (nM) v DMSO", y = "Log2FC (STC-15/DMSO) in run-on", 
       fill = "STC-15 (nM)")+
  theme_bw() + coord_cartesian(ylim = c(-1,2.5))

ggsave(plot = dose_response_intron_retention_boxplot, 
       filename = "analysis/rnaseq/plots/dose_response_stc15_caov3_intron_retention_boxplot.png", 
       width = 10, height = 10, units = "cm", dpi = 300) 

ggsave(plot = dose_response_runon_boxplot, 
       filename = "analysis/rnaseq/plots/dose_response_stc15_caov3_runon_boxplot.png", 
       width = 10, height = 10, units = "cm", dpi = 300) 

################################################################################
################################################################################
#### MC38 in vivo
################################################################################
#### IR

mc38_intron_retention_feature_id = mc38_invivo_intron_retention_results_table %>% 
  dplyr::filter(l2fc > 0.5, padj < 0.001, l2fc_de > 0, padj_de < 0.05, width > 250, ue_width >= 300) %>% pull(feature_id)
mc38_invivo_intron_retention_results_table_sigany = mc38_invivo_intron_retention_results_table %>% 
  filter(feature_id %in% mc38_intron_retention_feature_id) %>%
  mutate(time = factor(paste0(stc15_time_v_vehicle, "h"), levels = c("1h","4h","12h","24h","48h","120h"))) %>%
  mutate(cell = ifelse(grepl("tumour", cell), "Tumour", "Blood"))

data_table_mc38_invivo_intron_retention_results_table_sigany = mc38_invivo_intron_retention_results_table_sigany
names(data_table_mc38_invivo_intron_retention_results_table_sigany)[c(10:23)] = c("upstream_exon_range", "upstream_exon_type", "upstream_exon_width",
                                                                                  "downstream_exon_range", "downstream_exon_type", "downstream_exon_width",
                                                                                  "exon_base_mean", "DEXseq_padj", "DEXseq_log2FC", "DESeq_l2fc", "DESeq_padj", 
                                                                                  "stc15_time_v_vehicle", "tissue", "time")

write.csv(data_table_mc38_invivo_intron_retention_results_table_sigany, 
          file = "analysis/rnaseq/tables/data_table_significant_intron_retention_events_in_blood_or_tumour_at_any_timepoint_300nt_upstream_exon_widths.csv", quote = F, row.names = F)

mc38_invivo_intron_retention_results_table_sigany_boxplot = 
  ggplot(mc38_invivo_intron_retention_results_table_sigany, aes(x = time, y = l2fc, fill = cell))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.2, size = 0.5)+
  facet_grid(~cell, scales = "free_x")+ 
  theme_bw()+
  labs(x = "Time (h) after STC-15 dose v vehicle", 
       y = "Log2 FC (Time after last STC-15 dose/vehicle)\nin intron retention (significant events for any comparison)",
       fill = "Tissue type")+
  coord_cartesian(ylim = c(-0.5, 1.5))

ggsave(mc38_invivo_intron_retention_results_table_sigany_boxplot, 
       filename = "analysis/rnaseq/plots/mc38_invivo_intron_retention_results_table_sigany_boxplot.png",
       height = 12, width = 17, dpi = 300, units = "cm")

################################################################################
#### RO

mc38_runon_gene_id = mc38_invivo_runon_results_table %>% 
  dplyr::filter(l2fc > 1, padj < 0.00001, gene_l2fc < 1) %>% pull(gene_id) %>% unique()
mc38_invivo_runon_results_table_sigany = mc38_invivo_runon_results_table %>% 
  filter(gene_id %in% mc38_runon_gene_id) %>%
  mutate(time = factor(paste0(stc15_time_v_vehicle, "h"), levels = c("1h","4h","12h","24h","48h","120h"))) %>%
  mutate(cell = ifelse(grepl("tumour", cell), "Tumour", "Blood"))

data_table_mc38_invivo_runon_results_table_sigany = mc38_invivo_runon_results_table_sigany
names(data_table_mc38_invivo_runon_results_table_sigany)[3:15] = c("runon_chr", "runon_start", "runon_end", "runon_strand", "runon_width", "runon_base_mean", "runon_log2FC", "runon_padj", "gene_log2FC", "gene_padj", "stc15_time_v_vehicle", "tissue", "time")

write.csv(data_table_mc38_invivo_runon_results_table_sigany, file = "analysis/rnaseq/tables/data_table_significant_runon_events_in_blood_or_tumour_at_any_timepoint.csv", quote = F, row.names = F)

mc38_invivo_runon_results_table_sigany_boxplot = ggplot(mc38_invivo_runon_results_table_sigany, aes(x = time, y = l2fc, fill = cell))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(width = 0.1, alpha = 0.1, size = 0.5)+
  facet_grid(~cell, scales = "free_x")+ 
  theme_bw()+
  labs(x = "Time (h) after STC-15 dose v vehicle", 
       y = "Log2 FC (Time after last STC-15 dose/vehicle)\nin run-on (significant events for any comparison)",
       fill = "Tissue type")+
  coord_cartesian(ylim = c(-1, 4))

ggsave(mc38_invivo_runon_results_table_sigany_boxplot, 
       filename = "analysis/rnaseq/plots/mc38_invivo_runon_results_table_sigany_boxplot.png",
       height = 12, width = 17, dpi = 300, units = "cm")

################################################################################