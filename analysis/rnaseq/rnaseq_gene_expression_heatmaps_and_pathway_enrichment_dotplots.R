library(ComplexHeatmap)
library(circlize)
library(paletteer)
library(msigdbr)

human_neg_reg_viral_process_genes = msigdbr::msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
  filter(gs_name == "GOBP_NEGATIVE_REGULATION_OF_VIRAL_PROCESS") %>% pull(gene_symbol) %>% unique()
mouse_neg_reg_viral_process_genes = msigdbr::msigdbr(species = "Mus musculus", category = "C5", subcategory = "GO:BP") %>% 
  filter(gs_name == "GOBP_NEGATIVE_REGULATION_OF_VIRAL_PROCESS") %>% pull(gene_symbol)

################################################################################
####Log2FoldChange Heatmap plots

load("processing/rnaseq/multi_cell_line_stm3675/multi_cell_gene_expression_results_shrink_table.Rda")
load("processing/rnaseq/caov3_stc15_doseresponse/dose_response_gene_expression_results_shrink_table.Rda")
load("processing/rnaseq/mc38_invivo_model_stc15_timeseries/mc38_invivo_gene_expression_results_shrink_table.Rda")

write.csv(multi_cell_gene_expression_results_shrink_table, file = "analysis/rnaseq/tables/data_table_multi_cell_gene_expression_differential_expression_results.csv", quote = F, row.names = F)
write.csv(dose_response_gene_expression_results_shrink_table, file = "analysis/rnaseq/tables/data_table_dose_response_gene_expression_differential_expression_results.csv", quote = F, row.names = F)
write.csv(mc38_invivo_gene_expression_results_shrink_table, file = "analysis/rnaseq/tables/data_table_mc38_invivo_model_gene_expression_differential_expression_results.csv", quote = F, row.names = F)

################################################################################
#### Multi cell heatmap
################################################################################

cell_l2fc_table = multi_cell_gene_expression_results_shrink_table %>% 
  filter(gene_name %in% human_neg_reg_viral_process_genes) %>%
  mutate(l2fc = ifelse(padj < 1, log2FoldChange, 0)) %>%
  mutate(l2fc = ifelse(!is.na(l2fc), l2fc, 0)) %>%
  pivot_wider(id_cols = gene_name, names_from = cell, values_from = l2fc)

cell_changing_genes = cell_l2fc_table %>% column_to_rownames(var = "gene_name") %>% abs() %>% rowSums() != 0
cell_l2fc_matrix = cell_l2fc_table %>% filter(cell_changing_genes) %>% column_to_rownames(var = "gene_name") %>% as.matrix()

cell_l2fc_hm = Heatmap(cell_l2fc_matrix, column_labels = colnames(cell_l2fc_matrix),
                       row_labels = row.names(cell_l2fc_matrix), row_title = "GO:BP Negative regulation of viral process",
                       column_title = "Cells",
                       col = colorRamp2(seq(-2, 2, length.out = 30), rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))),
                       cluster_columns = T, cluster_rows = T, show_row_dend = F, show_column_dend = F,
                       column_names_side = "top", row_names_side = "left",
                       column_names_rot = 90,
                       heatmap_legend_param = list(title = "log2\n(STM3675/\nDMSO)\nFC",
                                                   legend_height = unit(4, "cm"), just = c("right", "top"), 
                                                   at = seq(-2,2,1), labels =  c("<-2","-1","0","1",">2")), 
                       border = "black", show_heatmap_legend = T, rect_gp = gpar(col = "grey60", lwd = 0.1))

png(filename = "analysis/rnaseq/plots/heatmap_cells_negative_regulation_viral_process.png",
    height = 34, width = 10, units = "cm", res = 300)
ComplexHeatmap::draw(cell_l2fc_hm)
################################################################################
dev.off()

################################################################################
#### Dose response Caov3 heatmap
################################################################################

dose_l2fc_table = dose_response_gene_expression_results_shrink_table %>% 
  filter(gene_name %in% human_neg_reg_viral_process_genes) %>%
  mutate(l2fc = ifelse(padj < 1, log2FoldChange, 0)) %>%
  mutate(l2fc = ifelse(!is.na(l2fc), l2fc, 0)) %>%
  pivot_wider(id_cols = gene_name, names_from = dose_v_0, values_from = l2fc) 

dose_changing_genes = dose_l2fc_table %>% column_to_rownames(var = "gene_name") %>% abs() %>% rowSums() != 0
dose_l2fc_matrix = dose_l2fc_table %>% filter(dose_changing_genes) %>% column_to_rownames(var = "gene_name") %>% as.matrix()

dose_l2fc_hm = Heatmap(dose_l2fc_matrix, column_labels = colnames(dose_l2fc_matrix),
                       row_labels = row.names(dose_l2fc_matrix), row_title = "GO:BP Negative regulation of viral process",
                       column_title = "STC-15 (nM)",
                       col = colorRamp2(seq(-2, 2, length.out = 30), rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))),
                       cluster_columns = F, cluster_rows = T, show_row_dend = F, show_column_dend = F,
                       column_names_side = "top", row_names_side = "left",
                       column_names_rot = 0, column_names_centered = T,
                       heatmap_legend_param = list(title = "log2\n(STC-15/\nDMSO)\nFC",
                                                   legend_height = unit(4, "cm"), just = c("right", "top"), 
                                                   at = seq(-2,2,1), labels =  c("<-2","-1","0","1",">2")), 
                       border = "black", show_heatmap_legend = T, rect_gp = gpar(col = "grey60", lwd = 0.1))

png(filename = "analysis/rnaseq/plots/heatmap_doses_negative_regulation_viral_process.png",
    height = 34, width = 11, units = "cm", res = 300)
ComplexHeatmap::draw(dose_l2fc_hm)
################################################################################
dev.off()

################################################################################
#### In vivo MC38 heatmap
################################################################################

mc38_invivo_gene_expression_results_shrink_table

tumour_l2fc_table = mc38_invivo_gene_expression_results_shrink_table %>% 
  filter(cell == "tumour") %>%
  filter(gene_name %in% mouse_neg_reg_viral_process_genes) %>%
  mutate(l2fc = ifelse(padj < 1, log2FoldChange, 0)) %>%
  mutate(l2fc = ifelse(!is.na(l2fc), l2fc, 0)) %>%
  pivot_wider(id_cols = gene_name, names_from = stc15_time_v_vehicle, values_from = l2fc)

blood_l2fc_table = mc38_invivo_gene_expression_results_shrink_table %>% 
  filter(cell == "blood") %>%
  filter(gene_name %in% mouse_neg_reg_viral_process_genes) %>%
  mutate(l2fc = ifelse(padj < 1, log2FoldChange, 0)) %>%
  mutate(l2fc = ifelse(!is.na(l2fc), l2fc, 0)) %>%
  pivot_wider(id_cols = gene_name, names_from = stc15_time_v_vehicle, values_from = l2fc)

################################################################################

changing_genes = tumour_l2fc_table %>% column_to_rownames(var = "gene_name") %>% abs() %>% rowSums() != 0 &
  blood_l2fc_table %>% column_to_rownames(var = "gene_name") %>% abs() %>% rowSums() != 0

tumour_l2fc_matrix = tumour_l2fc_table %>% filter(changing_genes) %>% column_to_rownames(var = "gene_name") %>% as.matrix()
blood_l2fc_matrix = blood_l2fc_table %>% filter(changing_genes) %>% column_to_rownames(var = "gene_name") %>% as.matrix()

ht = Heatmap(tumour_l2fc_matrix, cluster_rows = T, cluster_columns = F)
ht_drawn = draw(ht)
gene_order = row_order(ht_drawn)

tumour_l2fc_matrix_ordered = tumour_l2fc_matrix[gene_order,]
blood_l2fc_matrix_ordered = blood_l2fc_matrix[gene_order,]
combined_l2fc_matrix = cbind(blood_l2fc_matrix_ordered, tumour_l2fc_matrix_ordered)

################################################################################

invivo_l2fc_hm = Heatmap(combined_l2fc_matrix, column_labels = colnames(combined_l2fc_matrix),
                         column_split = c(rep("Time (h) after last dose\nBlood", 6), rep("Time (h) after last dose\nTumour", 6)),
                         row_labels = row.names(combined_l2fc_matrix), row_title = "GO:BP Negative regulation of viral process",
                         col = colorRamp2(seq(-2,2, length.out = 30), rev(paletteer_c("ggthemes::Red-Blue Diverging", 30))),
                         cluster_columns = F, cluster_rows = F, 
                         column_names_side = "top", row_names_side = "left",
                         column_names_rot = 0, column_names_centered = T, column_gap = unit(2, "mm"),
                         heatmap_legend_param = list(title = "log2\n(STC-15/\nVehicle)\nFC",
                                                     legend_height = unit(4, "cm"), just = c("right", "top"), 
                                                     at = seq(-2,2,1), labels =  c("<-2","-1","0","1",">2")), 
                         border = "black", show_heatmap_legend = T, rect_gp = gpar(col = "grey60", lwd = 0.1))


png(filename = "analysis/rnaseq/plots/heatmap_invivo_negative_regulation_viral_process.png",
    height = 30, width = 15, units = "cm", res = 300)
ComplexHeatmap::draw(invivo_l2fc_hm)
################################################################################
dev.off()

################################################################################
################################################################################
#### Pathway enrichment analysis

library(clusterProfiler)
library(org.Hs.eg.db)

cell_up_gene_list = multi_cell_gene_expression_results_shrink_table %>% filter(padj < 0.05, log2FoldChange > 0) %>%
  group_by(cell) %>% summarise(genes = list(gene_name), .groups = "drop") %>% deframe()
cell_down_gene_list = multi_cell_gene_expression_results_shrink_table %>% filter(padj < 0.05, log2FoldChange < 0) %>%
  group_by(cell) %>% summarise(genes = list(gene_name), .groups = "drop") %>% deframe()

dose_up_gene_list = dose_response_gene_expression_results_shrink_table %>% filter(padj < 0.05, log2FoldChange > 0) %>%
  group_by(dose_v_0) %>% summarise(genes = list(gene_name), .groups = "drop") %>% deframe()
dose_down_gene_list = dose_response_gene_expression_results_shrink_table %>% filter(padj < 0.05, log2FoldChange < 0) %>%
  group_by(dose_v_0) %>% summarise(genes = list(gene_name), .groups = "drop") %>% deframe()

cell_up_compare_result <- compareCluster(geneCluster = cell_up_gene_list,
                                         fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

cell_down_compare_result <- compareCluster(geneCluster = cell_down_gene_list,
                                           fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

dose_up_compare_result <- compareCluster(geneCluster = dose_up_gene_list,
                                         fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

dose_down_compare_result <- compareCluster(geneCluster = dose_down_gene_list,
                                           fun = "enrichGO", OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP")

options(scipen = 0)
ggsave(filename =  "analysis/rnaseq/plots/comparison_dotplot_cells_up_gene_expression.png",
       plot = dotplot(cell_up_compare_result, title = "", label_format = 50, showCategory = 5),
       height = 25, width = 25, units = "cm", dpi = 300)
ggsave(filename =  "analysis/rnaseq/plots/comparison_dotplot_cells_down_gene_expression.png",
       plot = dotplot(cell_down_compare_result, title = "", label_format = 50, showCategory = 5),
       height = 25, width = 25, units = "cm", dpi = 300)
ggsave(filename =  "analysis/rnaseq/plots/comparison_dotplot_doses_up_gene_expression.png",
       plot = dotplot(dose_up_compare_result, title = "", label_format = 50, showCategory = 5),
       height = 25, width = 25, units = "cm", dpi = 300)
ggsave(filename =  "analysis/rnaseq/plots/comparison_dotplot_doses_down_gene_expression.png",
       plot = dotplot(dose_down_compare_result, title = "", label_format = 50, showCategory = 5),
       height = 25, width = 25, units = "cm", dpi = 300)

################################################################################