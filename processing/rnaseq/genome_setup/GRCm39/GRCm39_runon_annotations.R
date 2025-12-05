library(rtracklayer)

################################################################################

gtf = rtracklayer::import("Mus_musculus.GRCm39.111.gtf")

################################################################################

ensembl_canonical_pct = gtf[gtf$tag %in% c("Ensembl_canonical") & gtf$gene_biotype == "protein_coding" & gtf$transcript_biotype == "protein_coding" &
                              gtf$type == "transcript" & !is.na(gtf$gene_name) & seqnames(gtf) %in% levels(seqnames(gtf))[1:21],] %>%
  as_tibble() %>% dplyr::group_by(gene_name) %>% filter(max(width) == width) %>% pull(transcript_id)

transcript = gtf[gtf$type == "transcript" & gtf$transcript_id %in% ensembl_canonical_pct]

redtxt = GenomicRanges::reduce(transcript)
fo1 = findOverlaps(redtxt, transcript)

redtxt$gene_id = as.vector(sapply(split(transcript[subjectHits(fo1)]$gene_id, queryHits(fo1)), 
                                  function(n) {n1 = 
                                    transcript[transcript$gene_id %in% n,]; 
                                  ifelse(strand(n1)[1] == "+", n1[which.max(end(n1)),]$gene_id, n1[which.min(start(n1)),]$gene_id)}))

rp = sort(redtxt[strand(redtxt) == "+"], decreasing = F)
rn = sort(redtxt[strand(redtxt) == "-"], decreasing = T)

################################################################################

runon_length = 20000

dsp = flank(rp, width= runon_length, start=F)
usp = flank(rp, width= 1, start=T)

dsn = flank(rn, width= runon_length, start=F)
usn = flank(rn, width= 1, start=T)

################################################################################

fop = findOverlaps(dsp, usp)
dsp1 = dsp
qh = queryHits(fop)
qhdups = duplicated(qh)
qh1 = qh[!qhdups]
sh = subjectHits(fop)
sh1 = sh[!qhdups]
end(dsp1[qh1]) = start(usp[sh1])

fon = findOverlaps(dsn, usn)
dsn1 = dsn
qh = queryHits(fon)
qhdups = duplicated(qh)
qh1 = qh[!qhdups]
sh = subjectHits(fon)
sh1 = sh[!qhdups]
start(dsn1[qh1]) = end(usn[sh1])

################################################################################
ro = sort(c(dsp1,dsn1))
ro$type = "exon"
start(ro)[start(ro) < 0] = 1
rtracklayer::export(ro, con = "Mmus_ensembl_canonical_protein_coding_runon.gtf")

