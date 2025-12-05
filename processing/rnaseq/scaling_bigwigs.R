################################################################################
#### Script for normalizing bigwig tracks

pos_bw_files = list.files("bigwig_directory/", pattern = "pos.bw", full.names = T)
neg_bw_files = gsub("pos", "neg", pos_bw_files)

bw_file_table = tibble(sample = gsub("[.]pos[.]bw", "", basename(pos_bw_files)), pos_bw_files, neg_bw_files)

"BigWig signal tracks were normalized by dividing each signal value by the total signal across the genome and multiplying by one billion. 
This yields a normalized signal representing signal-per-billion total signal units, facilitating cross-sample comparison."

for(a1 in 1:dim(bw_file_table)[1]){
  ##############################################################################
  sample_name = bw_file_table$sample[a1]
  bw_pos = bw_file_table$pos_bw_files[a1]
  bw_neg = bw_file_table$neg_bw_files[a1]
  print(sample_name)
  pos_gr = rtracklayer::import(bw_pos)
  neg_gr = rtracklayer::import(bw_neg)
  ##############################################################################
  total_score = sum(pos_gr$score) + sum(neg_gr$score)
  scaling_factor = 1e9/total_score
  pos_gr$score = round(pos_gr$score*scaling_factor, digits = 3)
  neg_gr$score = round(neg_gr$score*scaling_factor, digits = 3)
  rtracklayer::export(pos_gr, con = sprintf("bigwig_directory/scaled/%s_scaled.pos.bw", sample_name), format = "bigWig")
  rtracklayer::export(neg_gr, con = sprintf("bigwig_directory/scaled/%s_scaled.neg.bw", sample_name), format = "bigWig")
}

################################################################################
