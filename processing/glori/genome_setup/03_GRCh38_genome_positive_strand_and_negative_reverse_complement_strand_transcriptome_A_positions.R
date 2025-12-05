################################################################################
### Determine position of each A in protein coding and lncRNA transcriptome

library(GenomicRanges)
genome_fasta = "Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
gtf = rtracklayer::import("Homo_sapiens.GRCh38.108.gtf")

gf = Biostrings::readDNAStringSet(genome_fasta)
rcgf = reverseComplement(gf)
names1 = gsub(" .*", "", names(gf))

combined_genome = c(gf, rcgf)
names(combined_genome) = c(paste0("chrp_", names1), paste0("chrn_", names1))
cga = gregexpr("A", combined_genome)
cga1 = lapply(cga, function(n) as.vector(n))

combined_genome_names = names(combined_genome)
combined_genome_Alength = sapply(cga1, function(n) length(n))

cga1_chr = rep(combined_genome_names, times = combined_genome_Alength)
cga_df = data.frame(chr = cga1_chr, apos = unlist(cga1))
save(cga_df, file = "cga_df.Rda")

#load("cga_df.Rda")
annotations_gr = rtracklayer::import("GRCh38_annotations_combined_positive_and_negative_revcomp.gtf")
exons_gr = annotations_gr[annotations_gr$type == "exon",]
exons_gr1 = exons_gr[exons_gr$gene_biotype %in% c("protein_coding","lncRNA"),]
rexons_gr1 = reduce(exons_gr1)

chr_names = names(combined_genome)[!grepl("KI|GL", names(combined_genome))]

options(scipen = 999)
for(a1 in 1:length(chr_names)){
  print(chr_names[a1])
  apos = as.vector(cga[[a1]])
  chr_apos_gr = makeGRangesFromDataFrame(data.frame(chr = chr_names[a1], start = apos, end = apos))
  fo1 = findOverlaps(chr_apos_gr, rexons_gr1)
  write.table(data.frame(chr_apos_gr[queryHits(fo1)])[,1:2], 
              file = sprintf("apos_exons_%s.txt", chr_names[a1]),
              row.names = FALSE, col.names = FALSE, sep = "\t", quote = F)
}