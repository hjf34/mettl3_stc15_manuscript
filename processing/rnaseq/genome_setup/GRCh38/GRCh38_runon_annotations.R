library(rtracklayer)

################################################################################

gtf = rtracklayer::import("Homo_sapiens.GRCh38.108.gtf")
mane_transcript = gtf[gtf$type == "transcript" & gtf$tag %in% "MANE_Select" & seqnames(gtf) %in% levels(seqnames(gtf))[1:24]]

id2name = data.frame(mcols(mane_transcript)) %>% dplyr::select(c(gene_id, gene_name))
row.names(id2name) = id2name$gene_id

################################################################################

reduced_mane_transcript = GenomicRanges::reduce(mane_transcript)
fo1 = findOverlaps(reduced_mane_transcript, mane_transcript)

reduced_mane_transcript$combined_gene_id = sapply(split(mane_transcript[subjectHits(fo1)]$gene_id, queryHits(fo1)), function(n) paste(sort(unique(n)), collapse=";"))
reduced_mane_transcript$gene_id = as.vector(sapply(split(mane_transcript[subjectHits(fo1)]$gene_id, queryHits(fo1)), function(n) {n1 = mane_transcript[mane_transcript$gene_id %in% n,]; 
ifelse(strand(n1)[1] == "+", n1[which.max(end(n1)),]$gene_id, n1[which.min(start(n1)),]$gene_id)}))
reduced_mane_transcript$gene_name = id2name[reduced_mane_transcript$gene_id,]$gene_name

reduced_mane_transcript1 = reduced_mane_transcript
reduced_mane_transcript1$upstream_gene = NA
reduced_mane_transcript1$downstream_gene = NA
reduced_mane_transcript1[!is.na(follow(reduced_mane_transcript1)),]$upstream_gene = reduced_mane_transcript1[follow(reduced_mane_transcript1)[!is.na(follow(reduced_mane_transcript1))],]$gene_id
reduced_mane_transcript1[!is.na(precede(reduced_mane_transcript1)),]$downstream_gene = reduced_mane_transcript1[precede(reduced_mane_transcript1)[!is.na(precede(reduced_mane_transcript1))],]$gene_id

################################################################################
################################################################################

positive_transcripts = reduced_mane_transcript1[strand(reduced_mane_transcript1) == "+"]
negative_transcripts = reduced_mane_transcript1[strand(reduced_mane_transcript1) == "-"]
positive_transcripts$ds = precede(positive_transcripts)
negative_transcripts$ds = precede(negative_transcripts)

positive_transcripts$runon_final_limit = NA
positive_transcripts[!is.na(positive_transcripts$ds),]$runon_final_limit = start(positive_transcripts[positive_transcripts$ds[!is.na(positive_transcripts$ds)],]) - 1
negative_transcripts$runon_final_limit = NA
negative_transcripts[!is.na(negative_transcripts$ds),]$runon_final_limit = end(negative_transcripts[negative_transcripts$ds[!is.na(negative_transcripts$ds)],]) + 1

################################################################################

runon_length = 100000
tiling_length = 1000

positive_transcripts_runon_full_length = flank(positive_transcripts, width = runon_length, start = F)
end(positive_transcripts_runon_full_length) = 
  ifelse(end(positive_transcripts_runon_full_length) > positive_transcripts_runon_full_length$runon_final_limit & !is.na(positive_transcripts_runon_full_length$runon_final_limit), 
         positive_transcripts_runon_full_length$runon_final_limit, 
         end(positive_transcripts_runon_full_length))
positive_transcripts_runon_full_length1 = positive_transcripts_runon_full_length
end(positive_transcripts_runon_full_length1) = end(positive_transcripts_runon_full_length) - width(positive_transcripts_runon_full_length) %% tiling_length
positive_transcripts_runon_chunked = unlist(tile(positive_transcripts_runon_full_length1, width = tiling_length))

positive_transcripts_runon_chunked$gene_id = NA
fo_positive_transcripts_runon_chunked_hits = findOverlaps(positive_transcripts_runon_chunked, positive_transcripts_runon_full_length1)
positive_transcripts_runon_chunked$gene_id[queryHits(fo_positive_transcripts_runon_chunked_hits)] = positive_transcripts_runon_full_length1$gene_id[subjectHits(fo_positive_transcripts_runon_chunked_hits)]
positive_transcripts_runon_chunked$chunk_number = as.numeric(ave(positive_transcripts_runon_chunked$gene_id, positive_transcripts_runon_chunked$gene_id, FUN = seq_along))

negative_transcripts_runon_full_length = flank(negative_transcripts, width = runon_length, start = F)
start(negative_transcripts_runon_full_length) = ifelse(start(negative_transcripts_runon_full_length) < negative_transcripts_runon_full_length$runon_final_limit & !is.na(negative_transcripts_runon_full_length$runon_final_limit), 
                                                       negative_transcripts_runon_full_length$runon_final_limit, 
                                                       start(negative_transcripts_runon_full_length))
negative_transcripts_runon_full_length1 = negative_transcripts_runon_full_length
start(negative_transcripts_runon_full_length1) = start(negative_transcripts_runon_full_length) + width(negative_transcripts_runon_full_length) %% tiling_length
negative_transcripts_runon_chunked = unlist(tile(negative_transcripts_runon_full_length1, width = tiling_length))

negative_transcripts_runon_chunked$gene_id = NA
fo_negative_transcripts_runon_chunked_hits = findOverlaps(negative_transcripts_runon_chunked, negative_transcripts_runon_full_length1)
negative_transcripts_runon_chunked$gene_id[queryHits(fo_negative_transcripts_runon_chunked_hits)] = negative_transcripts_runon_full_length1$gene_id[subjectHits(fo_negative_transcripts_runon_chunked_hits)]
negative_transcripts_runon_chunked = sort(negative_transcripts_runon_chunked, decreasing = T)
negative_transcripts_runon_chunked$chunk_number = as.numeric(ave(negative_transcripts_runon_chunked$gene_id, negative_transcripts_runon_chunked$gene_id, FUN = seq_along))

################################################################################
################################################################################

runon_chunk = c(positive_transcripts_runon_chunked, negative_transcripts_runon_chunked[!start(negative_transcripts_runon_chunked) < 0])
runon_chunk$chunk = paste(runon_chunk$gene_id, sprintf("%03d", runon_chunk$chunk_number), sep = "_")
runon_chunk$type = "exon"

################################################################################
################################################################################
################################################################################

rtracklayer::export(runon_chunk, con = "Hsap_mane_runon_chunk1000.gtf")
runon_chunk1000_annotations = rtracklayer::import("Hsap_mane_runon_chunk1000.gtf")

df1 = data.frame(runon_chunk)
df1$score = 0
df1$start = start(runon_chunk) - 1

df1x = df1[,c("seqnames","start","end","chunk","score","strand")]

options(scipen = 1000)

write.table(df1x, file = "Hsap_mane_runon_chunk1000.bed", quote=F, sep="\t", row.names=F, col.names = F)

################################################################################
################################################################################
### Run featureCounts on Hsap_mane_runon_chunk1000.gtf

runon_chunk1000_table = as_tibble(read.table("featurecounts_mane_runon_chunk1000", skip=1, header = T))
names1 = gsub(".umidedup.bam", "", grep("umidedup.bam", names(runon_chunk1000_table), value = T))
names(runon_chunk1000_table) = c("chunk_id", "chr", "start", "end", "strand", "width", names1)

################################################################################

cts = runon_chunk1000 %>% dplyr::select(7:length(.)) %>% data.frame()
row.names(cts) = runon_chunk1000$chunk_id
sample_df = data.frame(row.names = colnames(cts), condition = factor(gsub("[[:digit:]]", "", colnames(cts))))

dds = DESeqDataSetFromMatrix(countData = cts, colData = sample_df, design = ~condition)
dds = estimateSizeFactors(dds)
norm_counts = counts(dds, normalized = T)

################################################################################
################################################################################
### Runon chunk annotations with treated normalized count greater than 10 in any repeat

roc = left_join(as_tibble(runon_chunk1000_annotations) %>% dplyr::select(gene_id, chunk, chunk_number), norm_counts %>% as_tibble(rownames = "chunk")) %>%
  left_join(id2name) %>% relocate(gene_name, .after = gene_id) %>% pivot_longer(cols = c(5:length(.))) %>% mutate(chunk_number = as.numeric(chunk_number))

roc1 = roc %>% filter(grepl("D1", name)) %>% dplyr::rename(c("D1" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell)) %>% 
  left_join(roc %>% filter(grepl("D2", name)) %>% dplyr::rename(c("D2" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(roc %>% filter(grepl("D3", name)) %>% dplyr::rename(c("D3" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(roc %>% filter(grepl("T1", name)) %>% dplyr::rename(c("T1" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(roc %>% filter(grepl("T2", name)) %>% dplyr::rename(c("T2" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell))) %>%
  left_join(roc %>% filter(grepl("T3", name)) %>% dplyr::rename(c("T3" = value, "cell" = name)) %>% mutate(cell = gsub("..$", "", cell)))

count_threshold = 10
runon_chunk_before_two_consecutive_chunks_below_threshold = roc1 %>% arrange(gene_id, cell, chunk_number) %>%
  mutate(across(c(D1,D2,D3,T1,T2,T3), function(n) n > count_threshold)) %>% 
  rowwise() %>% mutate(flag = any(c(T1,T2,T3))) %>% group_by(cell, gene_id) %>%
  mutate(lag_flag = lag(flag),
         is_false_pair = !flag & !lag_flag) %>%
  filter(chunk_number == ifelse(any(is_false_pair), max((which(is_false_pair)[1] - 1), 1), n()))

################################################################################
################################################################################
### At least 10 chunks long to be considered a runon i.e. 10 kb

runon_chunk_end_caov3 = runon_chunk_before_two_consecutive_chunks_below_threshold %>% ungroup() %>% filter(cell == "C" & chunk_number >= 10) %>%
  dplyr::select(c(gene_id, chunk_number))

runon_chunk_end_allcell = runon_chunk_before_two_consecutive_chunks_below_threshold %>% group_by(gene_id, gene_name) %>%
  summarise(across(chunk_number, max)) %>% filter(chunk_number >= 10) %>% 
  dplyr::select(c(gene_id, chunk_number)) %>% ungroup()

################################################################################
################################################################################

runon_chunk_caov3 = runon_chunk_end_caov3 %>% 
  tidyr::uncount(chunk_number) %>% group_by(gene_id) %>% 
  mutate(chunk_number = as.character(row_number())) %>% ungroup() %>%
  left_join(as_tibble(runon_chunk1000_annotations)) %>%
  group_by(gene_id) %>%
  summarise(seqnames = dplyr::first(seqnames),
            strand = dplyr::first(strand),
            start = min(start),
            end = max(end)) %>%
  left_join(id2name)

runon_chunk_caov3_gr = makeGRangesFromDataFrame(runon_chunk_caov3, keep.extra.columns = T)
runon_chunk_caov3_gr$type = "exon"

rtracklayer::export(runon_chunk_caov3_gr, con = "Hsap_mane_runon_joined_chunk1000_above_caov3_nc10.gtf")

df1 = data.frame(runon_chunk_caov3_gr)
df1$score = 0
df1$start = start(runon_chunk_caov3_gr) - 1

df1x = df1[,c("seqnames","start","end","gene_id","score","strand")]

options(scipen = 1000)

write.table(df1x, file = "Hsap_mane_runon_joined_chunk1000_above_caov3_nc10.bed", quote=F, sep="\t", row.names=F, col.names = F)

################################################################################

runon_chunk_allcell = runon_chunk_end_allcell %>% 
  tidyr::uncount(chunk_number) %>% group_by(gene_id) %>% 
  mutate(chunk_number = as.character(row_number())) %>% ungroup() %>%
  left_join(as_tibble(runon_chunk1000_annotations)) %>%
  group_by(gene_id) %>%
  summarise(seqnames = dplyr::first(seqnames),
            strand = dplyr::first(strand),
            start = min(start),
            end = max(end)) %>%
  left_join(id2name)

runon_chunk_allcell_gr = makeGRangesFromDataFrame(runon_chunk_allcell, keep.extra.columns = T)
runon_chunk_allcell_gr$type = "exon"

rtracklayer::export(runon_chunk_allcell_gr, con = "Hsap_mane_runon_joined_chunk1000_above_allcells_nc10.gtf")

df1 = data.frame(runon_chunk_allcell_gr)
df1$score = 0
df1$start = start(runon_chunk_allcell_gr) - 1

df1x = df1[,c("seqnames","start","end","gene_id","score","strand")]

options(scipen = 1000)

write.table(df1x, file = "Hsap_mane_runon_joined_chunk1000_above_allcells_nc10.bed", quote=F, sep="\t", row.names=F, col.names = F)

################################################################################