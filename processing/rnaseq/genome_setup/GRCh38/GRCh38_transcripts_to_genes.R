###Create gene map reference for salmon where each transcript is mapped to its corresponding gene
library(rtracklayer)

gtf = rtracklayer::import("Homo_sapiens.GRCh38.108.gtf")

transcripts = gtf[gtf$type == "transcript",]
tdf = as.data.frame(transcripts)
t_to_g = tdf[,c("transcript_id","gene_id")]

write.table(t_to_g, file = "processing/rnaseq/genome_setup/GRCh38/transcripts_to_genes.txt", sep= "\t", quote =F, col.names = F, row.names=F)