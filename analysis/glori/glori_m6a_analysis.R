################################################################################
################################################################################
library(MASS)  # for kde2d
library(RColorBrewer)

load("processing/glori/m6a_proportion_tables/m6a_proportion_table_filter_ag20_m6a25.Rda")

m6a_proportion_filter = 0.25
ag_count_filter = 20

m6a_proportion_table_filter_ag20_m6a25_both_replicates_caov3 = m6a_proportion_table_filter_ag20_m6a25 %>% filter(cell %in% c("CAOV3")) %>% 
  filter(D_m6a_1 >= m6a_proportion_filter & D_ag_1 >= ag_count_filter & D_m6a_2 >= m6a_proportion_filter & D_ag_2 >= ag_count_filter) %>% 
  dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>% distinct()

m6a_proportion_table_filter_ag20_m6a25_mean_caov3 = m6a_proportion_table_filter_ag20_m6a25 %>% filter(cell %in% c("CAOV3")) %>% 
  filter(D_m6a >= m6a_proportion_filter & D_ag >= ag_count_filter) %>% 
  dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>% distinct()

################################################################################

caov3_DTrepeat_boxplot = m6a_proportion_table_filter_ag20_m6a25_both_replicates_caov3 %>% 
  dplyr::select(c(xpos, D_m6a_1, D_m6a_2, T_m6a_1, T_m6a_2)) %>%
  pivot_longer(cols = c(D_m6a_1, D_m6a_2, T_m6a_1, T_m6a_2), names_to = "sample", values_to = "m6a") %>%
  dplyr::mutate(sample = gsub("D", "DMSO", gsub("T", "1.0 µM\nSTM3675", gsub("_m6a_", "\nrep ", sample)))) %>%
  dplyr::mutate(sample = factor(sample, levels = unique(sample))) %>%
  dplyr::mutate(treat = grepl("DMSO", sample)) %>%
  ggplot(aes(x = sample, y = m6a, fill = treat))+
  geom_violin(show.legend = F, alpha = 0.5)+
  scale_fill_manual(values = c("red", "blue"))+
  geom_boxplot(notch = T, coef = 50, width = 0.1, fill = "white")+
  theme_bw()+
  labs(x = "Caov-3 cells", 
       y = sprintf("m6A proportion (%s sites)", dim(m6a_proportion_table_filter_ag20_m6a25_both_replicates_caov3)[1]))

ggsave(plot = caov3_DTrepeat_boxplot,
       filename =  "analysis/glori/genome/plots/caov3_DTrepeat_boxplot.png",
       height = 12, width = 10, units = "cm", dpi = 300)

################################################################################

caov3_D1vD2 <- m6a_proportion_table_filter_ag20_m6a25_both_replicates_caov3 %>% dplyr::select(c("D_m6a_1", "D_m6a_2"))
caov3_DvT <- m6a_proportion_table_filter_ag20_m6a25_mean_caov3 %>% dplyr::select(c("D_m6a", "T_m6a"))

# Estimate 2D density
get_density <- function(x, y, ...) {
  dens <- kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

caov3_D1vD2$density <- get_density(caov3_D1vD2$D_m6a_1, caov3_D1vD2$D_m6a_2, n = 100)
caov3_DvT$density <- get_density(caov3_DvT$D_m6a, caov3_DvT$T_m6a, n = 100)

caov3_D1vD2_m6a_proportion_scatter_plot = ggplot(caov3_D1vD2, aes(D_m6a_1, D_m6a_2, color = density)) +
  geom_point(alpha = 0.08, size = 1, show.legend = F) +
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  viridis::scale_color_viridis(option = "plasma")+
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed")+
  theme_bw()+
  labs(title = sprintf("Caov-3 repeat correlation (%s sites)", dim(caov3_D1vD2)[1]),
       x = "DMSO m6A proportion repeat 1", y = "DMSO m6A proportion repeat 2")

caov3_DvT_mean_m6a_proportion_scatter_plot = ggplot(caov3_DvT, aes(D_m6a, T_m6a, color = density)) +
  geom_point(alpha = 0.08, size = 1, show.legend = F) +
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  viridis::scale_color_viridis(option = "plasma")+
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed")+
  theme_bw()+
  labs(title = sprintf("Caov-3 cell correlation (%s sites)", dim(caov3_DvT)[1]),
       x = "Mean DMSO m6A proportion", y = "Mean 1.0 µM STM3675 m6A proportion")

ggsave(plot = caov3_D1vD2_m6a_proportion_scatter_plot,
       filename =  "analysis/glori/genome/plots/caov3_D1vD2_m6a_proportion_scatter_plot.png",
       height = 12, width = 12, units = "cm", dpi = 300)

ggsave(plot = caov3_DvT_mean_m6a_proportion_scatter_plot,
       filename =  "analysis/glori/genome/plots/caov3_DvT_mean_m6a_proportion_scatter_plot.png",
       height = 12, width = 12, units = "cm", dpi = 300)

##Stats
dim(caov3_DvT)[1] # Number of sites
sum(caov3_DvT$D_m6a - caov3_DvT$T_m6a > 0)/dim(caov3_DvT)[1] #Proportion of sites reduced
sum(caov3_DvT$T_m6a < 0.25)/dim(caov3_DvT)[1] #Proportion of sites reduced below 25 %

################################################################################
################################################################################

m6a_proportion_table_filter_ag20_m6a25_mean_a549caov3 = m6a_proportion_table_filter_ag20_m6a25 %>% 
  filter(cell %in% c("CAOV3", "A549")) %>% 
  dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>% distinct()

m6a_proportion_table_filter_ag20_m6a25_mean_a549caov3 = m6a_proportion_table_filter_ag20_m6a25_mean_a549caov3 %>% 
  filter(xpos %in% 
           (m6a_proportion_table_filter_ag20_m6a25_mean_a549caov3 %>% 
              filter(D_m6a >= m6a_proportion_filter & D_ag >= ag_count_filter) %>% 
              group_by(xpos) %>% summarise(count = n()) %>% filter(count == 2) %>% pull(xpos)))

m6a_proportion_table_filter_ag20_m6a25_a549caov3_wide = 
  m6a_proportion_table_filter_ag20_m6a25_mean_a549caov3 %>% dplyr::select(xpos, cell, D_m6a) %>%
  pivot_wider(names_from = cell, values_from = D_m6a)

m6a_proportion_table_filter_ag20_m6a25_a549caov3_wide$density <- get_density(m6a_proportion_table_filter_ag20_m6a25_a549caov3_wide$A549, m6a_proportion_table_filter_ag20_m6a25_a549caov3_wide$CAOV3, n = 100)

caov3_v_a549_m6a_proportion_scatter_plot = ggplot(m6a_proportion_table_filter_ag20_m6a25_a549caov3_wide, aes(x = A549, y = CAOV3, col = density))+
  geom_point(alpha = 0.08, size = 1, show.legend = F) +
  coord_fixed(xlim=c(0,1), ylim=c(0,1))+
  viridis::scale_color_viridis(option = "plasma")+
  geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed")+
  theme_bw()+
  labs(x = "Mean A549 DMSO m6A proportion", y = "Mean CAOV3 DMSO m6A proportion",
       title = sprintf("A549 v Caov-3 (%s shared sites)", dim(m6a_proportion_table_filter_ag20_m6a25_a549caov3_wide)[1]))

ggsave(plot = caov3_v_a549_m6a_proportion_scatter_plot,
       filename =  "analysis/glori/genome/plots/caov3_v_a549_m6a_proportion_scatter_plot.png",
       height = 12, width = 12, units = "cm", dpi = 300)

################################################################################
################################################################################

shared_sites_ag20_m6a25 = m6a_proportion_table_filter_ag20_m6a25 %>% filter(D_m6a >= 0.25 & D_ag >= 20) %>% 
  group_by(xpos, gene_id) %>% summarise(count = n()) %>% filter(count == 6) %>% pull(xpos) %>% unique()
m6a_proportion_table_shared_sites = m6a_proportion_table_filter_ag20_m6a25 %>% filter(xpos %in% shared_sites_ag20_m6a25)
m6a_proportion_table_shared_sites_distinct = m6a_proportion_table_shared_sites %>% 
  dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>% distinct()

m6a_proportion_shared_sites_6cell_boxplot = ggplot(m6a_proportion_table_shared_sites_distinct, aes(x = cell, y = D_m6a, fill = cell))+
  geom_violin(fill = "white")+
  geom_boxplot(notch = T, coef = 50, width = 0.2)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  labs(x = "Cell", y = sprintf("m6A proportion (%s shared sites)", dim(m6a_proportion_table_shared_sites_distinct)[1]/6), fill = "Cell")+
  coord_cartesian(ylim = c(0, 1))

ggsave(filename = "analysis/glori/genome/plots/m6a_proportion_shared_sites_6cell_boxplot.png",
       plot = m6a_proportion_shared_sites_6cell_boxplot, width = 10, height = 12.5, units = "cm", dpi = 300)

