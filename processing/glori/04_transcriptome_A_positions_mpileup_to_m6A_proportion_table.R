################################################################################
library(tidyverse)
library(GenomicRanges)

seqinfo1 = seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
seqinfo1_df = as_tibble(data.frame(seqinfo1)[1:25,], rownames = "chr")

glori_directory = "."
bamfile_list = read.table(paste0(glori_directory, "alignment/genome/processing/all_alignGenome_gtoa_bamfiles.txt"))
samples = gsub("_align.*","", basename(bamfile_list$V1))

################################################################################
################################################################################

combined_chunks_mpileup_files = list.files(path = "mpileup/genome/", pattern = "combined", full.names = T)
chromosome_names = gsub("combined_|_chunks.tsv", "", basename(combined_chunks_mpileup_files))

################################################################################
################################################################################

arate_agtotal_calculator = function(sample, mpileup, base_names){
  base_count = lapply(strsplit(mpileup[,grep(sample, names(mpileup)), drop=T], ","), as.numeric)
  ag_count = t(mapply(function(X, Y) {names(X) = Y; X[c("A","G")]}, X = base_count, Y = base_names))
  ag_count[is.na(ag_count)] = 0
  colnames(ag_count) = c("A","G")
  agtotal = ag_count[,"A"] + ag_count[,"G"]
  arate = ag_count[,"A"]/agtotal
  arate[is.nan(arate)] = 0
  return(data.frame(arate, agtotal))
}

mpileup_processor = function(mpileup_file){
  ##############################################################################
  mpileup1 = tibble(read.table(mpileup_file))
  names(mpileup1) = c("chr", "pos", "ref", "alt", samples)
  chr_name = mpileup1$chr[1]
  chr_name1 = gsub("n_|p_|T", "", chr_name)
  strand = ifelse(grepl("n_", chr_name), "-", "+")
  chr_length = seqinfo1_df %>% filter(chr == chr_name1) %>% pull(seqlengths)
  if(strand == "-"){pos = chr_length - mpileup1$pos + 1}
  if(strand == "+"){pos = mpileup1$pos}
  base_names1 = strsplit(paste(mpileup1$ref, mpileup1$alt, sep = ","), ",")
  ##############################################################################
  arate_agtotal_list = lapply(samples, function(n) arate_agtotal_calculator(n, mpileup = mpileup1, base_names = base_names1))
  arate_df = tibble(data.frame(sapply(arate_agtotal_list, function(n) n[,1])))
  agtotal_df = tibble(data.frame(sapply(arate_agtotal_list, function(n) n[,2])))
  names(arate_df) = paste0("a_", samples)
  names(agtotal_df) = paste0("ag_", samples)
  ##############################################################################
  aag = tibble(chr = chr_name1, pos, strand, arate_df, agtotal_df) %>% arrange(pos)
  aag_file = paste0(glori_directory, "/ag_count/chromosomes/", 
                    sub("combined_", "aag_", sub("_chunks.tsv$", ".Rda", basename(mpileup_file))))
  save(aag, file = aag_file)
  print(aag)
}

for(a1 in 1:length(combined_chunks_mpileup_files)){
  print(combined_chunks_mpileup_files[a1])
  mpileup_processor(combined_chunks_mpileup_files[a1])
}

################################################################################
################################################################################

lf_aag = list.files(paste0(glori_directory, "ag_count/chromosomes/"))

main_chromosomes = gsub("[.]Rda|aag_", "", lf_aag)
main_chromosomes_positive = grep("p_", main_chromosomes, value = T)

for(a1 in 1:length(main_chromosomes_positive)){
  ##############################################################################
  chr_name_pos = main_chromosomes_positive[a1]
  print(chr_name_pos)
  chr_name_neg = sub("p", "n", chr_name_pos)
  aag_file_pos = paste0("aag_", chr_name_pos, ".Rda")
  aag_file_neg = paste0("aag_", chr_name_neg, ".Rda")
  load(paste0(glori_directory, "ag_count/chromosomes/", aag_file_pos))
  aag_pos = aag
  load(paste0(glori_directory, "ag_count/chromosomes/", aag_file_neg))
  aag_neg = aag
  ##############################################################################
  aag = bind_rows(aag_pos, aag_neg)
  amax = aag %>% dplyr::select(grep("a_", names(.))) %>% mutate(amax = do.call(pmax, (.))) %>% pull(amax)
  aagx = aag %>% filter(amax != 0)
  ##############################################################################
  ag_values = aagx %>% dplyr::select(c(1:3, grep("ag_", names(aag1)))) %>% pivot_longer(cols = grep("ag_", names(.)), names_to = "sample", values_to = "ag") %>% pull(ag)
  m6aag = aagx %>% dplyr::select(c(1:3, grep("a_", names(aag1)))) %>% pivot_longer(cols = grep("a_", names(.)), names_to = "sample", values_to = "m6a") %>%
    mutate(sample = gsub("a_", "", sample)) %>% mutate(ag = ag_values)
  print(m6aag)
  save(m6aag, file = paste0(glori_directory, "ag_count/m6aag_chromosomes/m6aag_", gsub("p_", "", chr_name_pos), ".Rda"))
  ##############################################################################
}

################################################################################
################################################################################

m6aag_chromosome_files = list.files(paste0(glori_directory, "/ag_count/m6aag_chromosomes/", full.names = T))

for(a1 in 1:length(m6aag_chromosome_files)){
  m6aag_chromosome_file = m6aag_chromosome_files[a1]
  chr_name = gsub("m6aag_|.Rda", "", basename(m6aag_chromosome_file))
  print(chr_name)
  load(m6aag_chromosome_file)
  print("chromosome loaded")
  ##############################################################################
  m6aag1 = m6aag %>% filter(grepl("1$", sample)) %>% dplyr::rename(c("m6a_1" = "m6a", "ag_1" = "ag")) %>%
    left_join(sample_info %>% dplyr::select(code, cell, treat), by = c("sample" = "code2")) %>% relocate(c(cell), .after = sample) %>%
    bind_cols(m6aag %>% filter(grepl("2$", sample)) %>% dplyr::rename(c("m6a_2" = "m6a", "ag_2" = "ag")) %>% dplyr::select(c(m6a_2, ag_2))) %>%
    relocate(m6a_2, .after = m6a_1) %>% mutate(xpos = paste(chr, pos, sep = "_")) %>%
    relocate(xpos, .after = pos) %>% dplyr::select(-sample)
  m6aag1$m6a = m6aag1 %>% dplyr::select(m6a_1, m6a_2) %>% rowMeans()
  m6aag1$ag = m6aag1 %>% dplyr::select(ag_1, ag_2) %>% rowMeans()
  m6aag1 = m6aag1 %>% mutate(m6a_x = pmax(m6a_1, m6a_2)) %>% mutate(ag_x = pmax(ag_1, ag_2)) %>%
    relocate(c(m6a, m6a_x), .after = m6a_2)
  print(m6aag1)
  ##############################################################################
  m6a_proportion_chromosome = m6aag1 %>% filter(treat == "DMSO") %>% dplyr::rename(c("D_m6a_1" = "m6a_1", "D_m6a_2" = "m6a_2", "D_m6a" = "m6a", "D_m6a_x" = "m6a_x", "D_ag_1" = "ag_1", "D_ag_2" = "ag_2", "D_ag" = "ag", "D_ag_x" = "ag_x")) %>%
     left_join(m6aag1 %>% filter(treat == "STM3675") %>% dplyr::select(c(chr, pos, strand, cell, m6a_1, m6a_2, m6a, m6a_x, ag_1, ag_2, ag, ag_x)) %>% 
                dplyr::rename(c("T_m6a_1" = "m6a_1", "T_m6a_2" = "m6a_2", "T_m6a" = "m6a", "T_m6a_x" = "m6a_x", "T_ag_1" = "ag_1", "T_ag_2" = "ag_2", "T_ag" = "ag", "T_ag_x" = "ag_x"))) %>%
    dplyr::select(-treat)
  save(m6a_proportion_chromosome, file = paste0(glori_directory, "/ag_count/m6a_proportion_chromosomes/m6a_proportion_chromosome_", chr_name, ".Rda"))
  ##############################################################################
  m6a_proportion_chromosome_filtered = m6a_proportion_chromosome %>% filter(xpos %in% (m6a_proportion_chromosome %>% filter((D_m6a_1 >= 0.20 & D_ag_1 >= 10) | (D_m6a_2 >= 0.20 & D_ag_2 >= 10)) %>% pull(xpos) %>% unique()))
  print(m6a_proportion_chromosome_filtered)
  save(m6a_proportion_chromosome_filtered, file = paste0(glori_directory, "/ag_count/m6a_proportion_chromosomes/m6a_proportion_chromosome_filtered_m6a_p2_ag_10_", chr_name, ".Rda"))
}

################################################################################
################################################################################

m6a_proportion_chromosome_filtered_files = list.files(paste0(glori_directory, "/ag_count/m6a_proportion_chromosomes/"), full.names = T, pattern = "filtered")

m6a_proportion_list = list()
for(a1 in 1:length(m6a_proportion_chromosome_filtered_files)){
  print(a1)
  m6a_proportion_chromosome_filtered_file = m6a_proportion_chromosome_filtered_files[a1]
  load(m6a_proportion_chromosome_filtered_file)
  m6a_proportion_list[[a1]] = m6a_proportion_chromosome_filtered
}

m6a_proportion_table = do.call(rbind, m6a_proportion_list)
m6a_proportion_table$chr = factor(m6a_proportion_table$chr, levels = seqinfo1_df$chr[1:24])
m6a_proportion_table$strand = factor(m6a_proportion_table$strand, levels = c("+", "-"))

m6a_proportion_table = m6a_proportion_table %>% arrange(chr, strand, pos)

################################################################################

library(Biostrings)

dna <- readDNAStringSet("Homo_sapiens.GRCh38.dna.primary_assembly.fa")
gtf = rtracklayer::import("Homo_sapiens.GRCh38.108.gtf")

dna1 = dna[c(1:22,24:25)]
names(dna1) = paste0("chr", gsub(" .*", "", names(dna)[c(1:22,24:25)]))

m6a_proportion_table_motif = m6a_proportion_table %>% mutate(start = pos-2, end = pos+2) %>% relocate(c(start, end), .after = strand)
m6a_proportion_table_motif_5mer_gr = makeGRangesFromDataFrame(m6a_proportion_table_motif, start.field = "start", end.field = "end")
m6a_proportion_table_motif = m6a_proportion_table_motif %>% mutate(dna5mer = as.character(BSgenome::getSeq(dna1, m6a_proportion_table_motif_5mer_gr))) %>% 
  relocate(c(dna5mer), .after = end) %>% dplyr::select(-c(start, end))

################################################################################
################################################################################

pos_gr = m6a_proportion_table_motif %>% dplyr::select(c(1:4)) %>% makeGRangesFromDataFrame(start.field = "pos", end.field = "pos", keep.extra.columns = T)

mane_exons_gtf = gtf[gtf$type == "exon" & gtf$tag %in% "MANE_Select",c("gene_id","gene_name","exon_number","transcript_name")]

################################################################################

## Determine if exon is first, internal, last or only exon in gene (MANE_Select transcripts)
megt = as_tibble(mane_exons_gtf)
megt = megt %>% mutate(seqnames = paste0("chr", seqnames))
megt$exon_number = as.integer(megt$exon_number)
megt = megt %>% group_by(gene_id) %>% mutate(exon_type = ifelse(exon_number == max(exon_number), ifelse(exon_number == 1, "solo", "last"), ifelse(exon_number == 1, "first", "internal")))
megt$genomic_data = paste0(megt$seqnames, ":", megt$start, "-", megt$end, ":", megt$strand)
megt_gr = makeGRangesFromDataFrame(megt, keep.extra.columns = T)

################################################################################
################################################################################

fo1 = findOverlaps(pos_gr, megt_gr)
m6a_proportion_table_motif_with_mane_exons = bind_cols(m6a_proportion_table_motif[queryHits(fo1),], megt[subjectHits(fo1),c(6:8,10,11,4)]) %>% relocate(c(gene_id, gene_name, exon_number, exon_type, width, genomic_data), .after = strand)

m6a_proportion_table_filter_ag10_m6a20 = full_join(m6a_proportion_table_motif_with_mane_exons, m6a_proportion_table_motif)
m6a_proportion_table_filter_ag10_m6a20 = m6a_proportion_table_filter_ag10_m6a20 %>% arrange(chr, strand, pos) %>% 
  dplyr::rename(c("exon_annotation" = "genomic_data", "exon_width" = "width")) %>%
  relocate(exon_annotation, .after = gene_name)

any_sites_ag20_m6a25 = m6a_proportion_table_filter_ag10_m6a20 %>% filter((D_m6a_1 >= 0.25 & D_ag_1 >= 20) | (D_m6a_2 >= 0.25 & D_ag_2 >= 20)) %>% pull(xpos) %>% unique()
m6a_proportion_table_filter_ag20_m6a25 = m6a_proportion_table_filter_ag10_m6a20 %>% filter(xpos %in% any_sites_ag20_m6a25)
save(m6a_proportion_table_filter_ag20_m6a25, file = "processing/glori/m6a_proportion_tables/m6a_proportion_table_filter_ag20_m6a25.Rda")

options(scipen = 999)
data_table_m6a_proportion_table_filter_ag20_m6a25 = m6a_proportion_table_filter_ag20_m6a25
names(data_table_m6a_proportion_table_filter_ag20_m6a25)[13:28] = c("DMSO_m6a_rep1", "DMSO_m6a_rep2", "DMSO_m6a_mean", "DMSO_m6a_max", "DMSO_AG_rep1", "DMSO_AG_rep2", "DMSO_AG_mean", "DMSO_AG_max",
                                                                    "STM3675_m6a_rep1", "STM3675_m6a_rep2", "STM3675_m6a_mean", "STM3675_m6a_max", "STM3675_AG_rep1", "STM3675_AG_rep2", "STM3675_AG_mean", "STM3675_AG_max")
write.csv(data_table_m6a_proportion_table_filter_ag20_m6a25, file = "processing/glori/m6a_proportion_tables/data_table_m6a_proportion_table_filter_ag20_m6a25.csv", quote = F, row.names = F)
zip(zipfile = "processing/glori/m6a_proportion_tables/data_table_m6a_proportion_table_filter_ag20_m6a25.csv.zip", files = "processing/glori/m6a_proportion_tables/data_table_m6a_proportion_table_filter_ag20_m6a25.csv")

################################################################################
