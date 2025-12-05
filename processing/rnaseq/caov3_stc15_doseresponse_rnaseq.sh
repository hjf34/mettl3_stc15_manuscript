### in vitro RNA-seq sample processing

wd="dose_response"

sample_info="processing/rnaseq/caov3_stc15_doseresponse/sample_info/doseresponse_sample_info.csv"
samples=$(cat $sample_info | cut -f1 -d, | tail -n+2 | uniq)

raw_fastq_directory="${wd}/raw_fastq/"
trimmed_fastq_directory="${wd}/trimmed_fastq/"
salmon_counts_directory="${wd}/salmon_counts/"
alignment_directory="${wd}/alignment"
alignment_star_directory="${alignment_directory}/STAR_bam"
bw_directory="${alignment_directory}/STAR_bw"

mkdir -p ${combined_fastq_directory}
mkdir -p ${trimmed_fastq_directory}
mkdir -p ${salmon_counts_directory}
mkdir -p ${alignment_directory}
mkdir -p ${alignment_star_directory}
mkdir -p ${bw_directory}

### Trimming fastq files

echo $samples
for sample in $samples
do
  echo $sample
  raw_fastq_file1="${raw_fastq_directory}/${sample}_1.fq.gz"
  raw_fastq_file2="${raw_fastq_directory}/${sample}_2.fq.gz"
  trim_galore -j 8 --paired ${combined_fastq_file1} ${combined_fastq_file2} -o ${trimmed_fastq_directory}
done

### Salmon counts

SalmonIndexDirectory="GRCh38_salmon_index"

echo $samples
for sample in $samples
do
  echo $sample
  echo "salmon transcript counting..."
  trimmed_fastq_file1="${trimmed_fastq_directory}/${sample}_1_val_1.fq.gz"
  trimmed_fastq_file2="${trimmed_fastq_directory}/${sample}_2_val_2.fq.gz"
  salmon_counts="${salmon_counts_directory}/${sample}"
  salmon quant -i $SalmonIndexDirectory -l A -1 ${trimmed_fastq_file1} -2 ${trimmed_fastq_file2} -p 16 -o $salmon_counts
done

### STAR align

STARIndexDirectory="GRCh38_STAR_index"

echo $samples
for sample in $samples
do
  echo $sample
  echo "STAR alignment..."
  trimmed_fastq_file1="${trimmed_fastq_directory}/${sample}_1_val_1.fq.gz"
  trimmed_fastq_file2="${trimmed_fastq_directory}/${sample}_2_val_2.fq.gz"
  align_star="${alignment_star_directory}/${sample}."
  STAR --runThreadN 16 --readFilesCommand zcat --genomeDir ${STARIndexDirectory} --readFilesIn ${trimmed_fastq_file1} ${trimmed_fastq_file2} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI AS NM MD --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix ${align_star}
done

################################################################################
################################################################################
### bigwig generation - strands are wrong (pos is neg and vice versa)

chromsizes="GRCh38_chromSizes"

echo $samples
for sample in $samples
do
  ##############################################################################
  echo $sample
  bam_file="${alignment_star_directory}/${sample}.Aligned.sortedByCoord.out.bam"
  samtools index ${bam_file}
  ls $bam_file
  bg_pos="${bw_directory}/${sample}.pos.bedGraph"
  bg_neg="${bw_directory}/${sample}.neg.bedGraph"
  bg_clip_pos="${bw_directory}/${sample}.pos_clip.bedGraph"
  bg_clip_neg="${bw_directory}/${sample}.neg_clip.bedGraph"
  bg_sort_pos="${bw_directory}/${sample}.pos_sort.bedGraph"
  bg_sort_neg="${bw_directory}/${sample}.neg_sort.bedGraph"
  bw_pos="${bw_directory}/${sample}.pos.bw"
  bw_neg="${bw_directory}/${sample}.neg.bw"
  bedtools genomecov -ibam $bam_file -bg -strand - -split -du > $bg_pos &
  bedtools genomecov -ibam $bam_file -bg -strand + -split -du > $bg_neg
  wait
  bedClip $bg_pos $chromsizes $bg_clip_pos &
  bedClip $bg_neg $chromsizes $bg_clip_neg
  wait
  bedSort $bg_clip_pos $bg_sort_pos &
  bedSort $bg_clip_neg $bg_sort_neg
  wait
  bedGraphToBigWig $bg_sort_pos $chromsizes $bw_pos &
  bedGraphToBigWig $bg_sort_neg $chromsizes $bw_neg
  wait
  rm $bg_clip_pos $bg_clip_neg $bg_sort_pos $bg_sort_neg $bg_pos $bg_neg
done

################################################################################
###Intron GTF/RO GTF

ir_featurecounts_directory="${wd}/feature_counts/intron_retention/"
ro_featurecounts_directory="${wd}/feature_counts/runon/"
mkdir -p ${ir_featurecounts_directory}
mkdir -p ${ro_featurecounts_directory}

IntronAnnotationGTF="Hsap_mane_select_introns.saf"

bamfiles=$(ls ${wd}/alignment/STAR_bam/*.bam)

featurecounts_introns="${ir_featurecounts_directory}/featurecounts_introns_doseresponse"

##Map to introns SAF
featureCounts -p --countReadPairs -s 2 -T 16 -F "SAF" -f -O -M -a ${IntronAnnotationGTF} -o ${featurecounts_introns} ${bamfiles}

mane_runon_joined_chunk_above_caov3_nc10_AnnotationsGTF="Hsap_mane_runon_joined_chunk1000_above_caov3_nc10.gtf"
featurecounts_mane_runon_joined_chunk_above_caov3_nc10="${ro_featurecounts_directory}/featurecounts_mane_runon_joined_chunk_above_caov3_nc10_doseresponse"

featureCounts -p --countReadPairs -s 2 -T 16 -O -M -a ${mane_runon_joined_chunk_above_caov3_nc10_AnnotationsGTF} -o ${featurecounts_mane_runon_joined_chunk_above_caov3_nc10} ${bamfiles}
