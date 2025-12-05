################################################################################
load(file = "processing/glori/m6a_proportion_tables/m6a_proportion_spike_in_table.Rda")

m6a_proportion_spike_in_true_positive_positions = m6a_proportion_spike_in_table %>% filter(m6a) %>% 
  mutate(expected_Arate = rep(c(0.05,0.2,0.5,0.8,1.0), length(m6a_proportion_spike_in_table %>% pull(code) %>% unique()))) %>%
  mutate(expected_Arate_x = factor(expected_Arate))

caov3_true_positive_positions_plot = 
  ggplot(m6a_proportion_spike_in_true_positive_positions %>% filter(cell == "CAOV3"), 
         aes(expected_Arate, Arate, col = treatx))+
  geom_abline(lwd = 0.1, lty = 2)+
  geom_point(alpha = 0.5)+
  coord_fixed(xlim= c(0,1), ylim = c(0,1))+
  labs(x = "Expected m6A proportion", y = "Measured m6A proportion\n(A count) / (A + G count)", 
       col = "Treatment")+
  theme_bw()+
  scale_color_manual(values = c("blue", "red"))+
  scale_y_continuous(breaks=c(0.05,0.2,0.5,0.8,1.0),labels=c(0.05,0.2,0.5,0.8,1.0))+
  scale_x_continuous(breaks=c(0.05,0.2,0.5,0.8,1.0),labels=c(0.05,0.2,0.5,0.8,1.0))

multi_cell_true_positive_positions_plot = 
  ggplot(m6a_proportion_spike_in_true_positive_positions %>% filter(treat == "DMSO"), 
         aes(expected_Arate, Arate, col = cell))+
  geom_abline(lwd = 0.1, lty = 2)+
  geom_point(alpha = 0.5)+
  coord_fixed(xlim= c(0,1), ylim = c(0,1))+
  labs(x = "Expected m6A proportion", y = "Measured m6A proportion\n(A count) / (A + G count)", 
       col = "Cell")+
  theme_bw()+
  scale_y_continuous(breaks=c(0.05,0.2,0.5,0.8,1.0),labels=c(0.05,0.2,0.5,0.8,1.0))+
  scale_x_continuous(breaks=c(0.05,0.2,0.5,0.8,1.0),labels=c(0.05,0.2,0.5,0.8,1.0))

ggsave("analysis/glori/spike_in/plots/caov3_true_positive_positions_plot.png", 
       plot = caov3_true_positive_positions_plot, 
       height = 12, width = 12, units = "cm", dpi = 300)
ggsave("analysis/glori/spike_in/plots/multi_cell_true_positive_positions_plot.png", 
       plot = multi_cell_true_positive_positions_plot, 
       height = 12, width = 12, units = "cm", dpi = 300)

################################################################################