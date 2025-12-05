################################################################################

seqinfo1 = seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
seqinfo1_df = as_tibble(data.frame(seqinfo1)[1:25,], rownames = "chr")

bigwig_directory = "processing/glori/bigwigs/"

load(file = "processing/glori/m6a_proportion_tables/m6a_proportion_table_filter_ag20_m6a25.Rda")

m6at = m6a_proportion_table_filter_ag20_m6a25 %>% dplyr::select(-c(gene_id, gene_name, exon_annotation, exon_number, exon_type, exon_width, dna5mer)) %>%
  distinct()
cells = m6at %>% pull(cell) %>% unique()

################################################################################
#### Untreated bigwigs

for(a1 in 1:length(cells)){
  cell1 = cells[a1]
  ##############################################################################
  m6at_c = m6at %>% filter(cell == cell1, D_m6a >= 0.25, D_ag >= 20) %>% 
    dplyr::rename("score" = D_m6a) %>% dplyr::select(chr, pos, strand, score) %>% mutate(score = round(score, digits = 2)) %>%
    makeGRangesFromDataFrame(start.field = "pos", end.field = "pos", keep.extra.columns = T, seqinfo = seqinfo1)
  m6at_cp = m6at_cb[strand(m6at_c) == "+"]
  m6at_cn = m6at_cb[strand(m6at_c) == "-"]
  file_name_stem = paste0(bigwig_directory, sprintf("mean_DMSO_m6A_proportion_ag20_m6ap25_%s", cell1))
  file_pos = paste0(file_name_stem, "_pos.bw")
  file_neg = paste0(file_name_stem, "_neg.bw")
  rtracklayer::export.bw(m6at_cp, con=file_pos)
  rtracklayer::export.bw(m6at_cn, con=file_neg)
  ##############################################################################
}

################################################################################
################################################################################
### Treated bigwigs

cell1 = "CAOV3"
##############################################################################
m6at_c = m6at %>% filter(cell == cell1, D_m6a >= 0.25, D_ag >= 20) %>% 
  dplyr::rename("score" = T_m6a) %>% dplyr::select(chr, pos, strand, score) %>% mutate(score = round(score, digits = 2)) %>%
  makeGRangesFromDataFrame(start.field = "pos", end.field = "pos", keep.extra.columns = T, seqinfo = seqinfo1)
m6at_cp = m6at_cb[strand(m6at_c) == "+"]
m6at_cn = m6at_cb[strand(m6at_c) == "-"]
file_name_stem = paste0(bigwig_directory, sprintf("mean_STM3675_m6A_proportion_ag20_m6ap25_%s", cell1))
file_pos = paste0(file_name_stem, "_pos.bw")
file_neg = paste0(file_name_stem, "_neg.bw")
rtracklayer::export.bw(m6at_cp, con=file_pos)
rtracklayer::export.bw(m6at_cn, con=file_neg)
##############################################################################
