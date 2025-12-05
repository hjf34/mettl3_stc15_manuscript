library(rtracklayer)
library(tidyverse)

################################################################################

gtf = rtracklayer::import("Mus_musculus.GRCm39.111.gtf")

library(GenomicRanges)
ensembl_canonical_pct = gtf[gtf$tag %in% c("Ensembl_canonical") & gtf$gene_biotype == "protein_coding" & gtf$transcript_biotype == "protein_coding" &
                              gtf$type == "transcript" & !is.na(gtf$gene_name) & seqnames(gtf) %in% levels(seqnames(gtf))[1:21],] %>%
  as_tibble() %>% dplyr::group_by(gene_name) %>% filter(max(width) == width) %>% pull(transcript_id)

gr = gtf[gtf$type == "exon" & gtf$transcript_id %in% ensembl_canonical_pct]

# find gaps (introns) between exons split by gene
s1 = split(gr, gr$gene_id)
gaps1 = lapply(s1, function(n) gaps(n, min(start(n))))
L1 = sapply(gaps1, length)
gaps2 = gaps1[L1 != 0] #exclude intronless genes

names_gaps2 = names(gaps2)
gaps3 = gaps2
for(a1 in 1:length(gaps2)){
  g1 = gaps2[[a1]]
  g1$gene_id = names_gaps2[a1]
  gaps3[[a1]] = g1
}

gaps4 = lapply(gaps3, as.data.frame)
grx = do.call(rbind, gaps4)

row.names(grx) = NULL
introns = grx[,c(6,1,2,3,5)]
names(introns) = c("GeneID", "Chr", "Start", "End", "Strand")

write.table(introns, file="Mmus_ensembl_canonical_protein_coding_introns.saf", quote = F, sep="\t", row.names = F )

################################################################################
#### Annotate intron annotations with type and width of upstream exon

introns = as_tibble(read.table(file="Mmus_ensembl_canonical_protein_coding_introns.saf", sep="\t", header = T)) 
introns1 = introns %>%
  arrange(GeneID, Start) %>% group_by(GeneID) %>% mutate(intron = paste0("e", seq(n())))
introns_gr = makeGRangesFromDataFrame(introns1, keep.extra.columns = T)

ensembl_canonical_pct = gtf[gtf$tag %in% c("Ensembl_canonical") & gtf$gene_biotype == "protein_coding" & gtf$transcript_biotype == "protein_coding" &
                              gtf$type == "transcript" & !is.na(gtf$gene_name) & seqnames(gtf) %in% levels(seqnames(gtf))[1:21],] %>%
  as_tibble() %>% dplyr::group_by(gene_name) %>% filter(max(width) == width) %>% pull(transcript_id)

exons_gr = gtf[gtf$type == "exon" & gtf$transcript_id %in% ensembl_canonical_pct]

egt = as_tibble(exons_gr)
egt$exon_number = as.integer(egt$exon_number)
egt = egt %>% group_by(gene_id) %>% mutate(exon_type = ifelse(exon_number == max(exon_number), ifelse(exon_number == 1, "solo", "last"), ifelse(exon_number == 1, "first", "internal")))
egt$gdata = paste0("chr", egt$seqnames, ":", egt$start, "-", egt$end, ":", egt$strand)
egt_gr = makeGRangesFromDataFrame(egt, keep.extra.columns = T)

################################################################################
################################################################################

follow1 = follow(introns_gr, egt_gr)
introns1$upstream_exon_gene = egt_gr[follow1,]$gene_id
introns1$upstream_exon_type = egt_gr[follow1,]$exon_type
introns1$upstream_exon_width = as_tibble(egt_gr[follow1,0])$width
introns1$upstream_range = as_tibble(egt_gr[follow1,0]) %>% mutate(gdata = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>% pull(gdata)

introns1_noerror = introns1[introns1$GeneID == introns1$upstream_exon_gene,]
introns1_error = introns1[introns1$GeneID != introns1$upstream_exon_gene,]

for(a1 in 1:dim(introns1_error)[1]){
  img1 = makeGRangesFromDataFrame(introns1_error[a1,])
  egtx = egt_gr[egt_gr$gene_id == introns1_error[a1,]$GeneID]
  egt1 = egtx[follow(img1, egtx),]
  introns1_error[a1,]$upstream_exon_gene = egt1$gene_id
  introns1_error[a1,]$upstream_exon_type = egt1$exon_type
  introns1_error[a1,]$upstream_exon_width = as_tibble(egt1[,0])$width
  introns1_error[a1,]$upstream_range = as_tibble(egt1[,0]) %>% mutate(gdata = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>% pull(gdata)
}

introns1 = bind_rows(introns1_error, introns1_noerror) %>% ungroup() %>%
  arrange(GeneID, Start)

################################################################################

precede1 = precede(introns_gr, egt_gr)
introns1$downstream_exon_gene = egt_gr[precede1,]$gene_id
introns1$downstream_exon_type = egt_gr[precede1,]$exon_type
introns1$downstream_exon_width = as_tibble(egt_gr[precede1,0])$width
introns1$downstream_range = as_tibble(egt_gr[precede1,0]) %>% mutate(gdata = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>% pull(gdata)

introns1_noerror = introns1[introns1$GeneID == introns1$downstream_exon_gene,]
introns1_error = introns1[introns1$GeneID != introns1$downstream_exon_gene,]

for(a1 in 1:dim(introns1_error)[1]){
  img1 = makeGRangesFromDataFrame(introns1_error[a1,])
  egtx = egt_gr[egt_gr$gene_id == introns1_error[a1,]$GeneID]
  egt1 = egtx[precede(img1, egtx),]
  introns1_error[a1,]$downstream_exon_gene = egt1$gene_id
  introns1_error[a1,]$downstream_exon_type = egt1$exon_type
  introns1_error[a1,]$downstream_exon_width = as_tibble(egt1[,0])$width
  introns1_error[a1,]$downstream_range = as_tibble(egt1[,0]) %>% mutate(gdata = paste0("chr", seqnames, ":", start, "-", end, ":", strand)) %>% pull(gdata)
}

introns1 = bind_rows(introns1_error, introns1_noerror) %>% ungroup() %>%
  arrange(GeneID, Start)

################################################################################

introns1$featureID = factor(introns1$intron, levels = unique(introns1$intron))
intron_with_updn_exon = introns1 %>% arrange(GeneID, featureID)

names(intron_with_updn_exon) = c("gene_id", "chr", "start", "end", "strand", "intron",
                                 "ue_gene_id", "ue_type", "ue_width", "ue_range", 
                                 "de_gene_id", "de_type", "de_width", "de_range", "feature_id")
intron_with_updn_exon = intron_with_updn_exon %>% dplyr::select(c(gene_id, feature_id, chr, start, end, strand, ue_range, ue_type, ue_width, de_range, de_type, de_width))
intron_with_updn_exon = intron_with_updn_exon %>% left_join(egt %>% dplyr::select(c(gene_id, gene_name)) %>% distinct()) %>% relocate(gene_name, .after = gene_id)
mmus_ensembl_canonical_protein_coding_introns_with_updn_exons = intron_with_updn_exon

save(mmus_ensembl_canonical_protein_coding_introns_with_updn_exons, file = "mmus_ensembl_canonical_protein_coding_introns_with_updn_exons.Rda")
