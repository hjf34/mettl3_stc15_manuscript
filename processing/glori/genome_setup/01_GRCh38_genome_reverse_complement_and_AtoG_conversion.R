################################################################################
library(Biostrings)
library(tidyr)
library(GenomicRanges)
library(rtracklayer)

# Download GRCh38 genome and annotation file from Ensembl https://ftp.ensembl.org/
genome_fasta = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf = rtracklayer::import("Homo_sapiens.GRCh38.108.gtf")

gf = Biostrings::readDNAStringSet(genome_fasta)
rcgf = reverseComplement(gf)

gf_AG = DNAStringSet(gsub("A", "G", gf))
rcgf_AG = DNAStringSet(gsub("A", "G", rcgf))

################################################################################

names1 = gsub(" .*", "", names(gf))

combined_genome_AG = c(gf_AG, rcgf_AG)
names(combined_genome_AG) = c(paste0("chrp_", names1), paste0("chrn_", names1))
grch38_combined_strands_AtoG_fasta = "GRCh38_DNA_combined_strands_AtoG.fa"
writeXStringSet(combined_genome_AG, grch38_combined_strands_AtoG_fasta, compress = F, width = 60)

################################################################################

positive_annotations = gtf[strand(gtf) == "+",]
negative_annotations = gtf[strand(gtf) == "-",]

chromosome_lengths = data.frame(name = gsub(" .*","", names(gf)), len = width(gf))

seqlevels1 = seqlevels(negative_annotations)
negative_annotations_list = list()
for(a1 in 1:length(seqlevels1)){
  sl1 = seqlevels1[a1]
  print(sl1)
  negative_annotations1 = negative_annotations[seqnames(negative_annotations) == sl1,]
  chromosome_lengths1 = chromosome_lengths[chromosome_lengths$name == sl1,]$len
  startx = chromosome_lengths1 - start(negative_annotations1) + 1
  endx = chromosome_lengths1 - end(negative_annotations1) + 1
  negative_annotations1_df = data.frame(negative_annotations1)
  negative_annotations1_df[,"start"] = endx
  negative_annotations1_df[,"end"] = startx
  negative_annotations_list[[a1]] = negative_annotations1_df
}

positive_annotations_df = as_tibble(positive_annotations) %>%
  mutate(gene_id = factor(gene_id, levels = unique(gene_id)), transcript_id = factor(transcript_id, levels = unique(transcript_id)),
         seqnames = paste0("chrp_", as.vector(seqnames))) %>% 
  dplyr::arrange(seqnames, gene_id, desc(is.na(transcript_id)), transcript_id, start, type)
  
negative_annotations_df = as_tibble(do.call(rbind, negative_annotations_list)) %>% mutate(strand = "+") %>%
  mutate(gene_id = factor(gene_id, levels = unique(gene_id)), transcript_id = factor(transcript_id, levels = unique(transcript_id)),
         seqnames = paste0("chrn_", as.vector(seqnames)))  %>% 
  dplyr::arrange(seqnames, gene_id, desc(is.na(transcript_id)), transcript_id, start, type)

################################################################################
################################################################################

annotations_gr =  makeGRangesFromDataFrame(rbind(positive_annotations_df, negative_annotations_df), keep.extra.columns = T)
rtracklayer::export(object = annotations_gr, con = "GRCh38_annotations_combined_positive_and_negative_revcomp.gtf")

################################################################################
