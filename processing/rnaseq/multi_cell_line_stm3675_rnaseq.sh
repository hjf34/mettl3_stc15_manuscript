### in vitro RNA-seq sample processing

STARIndexDirectory="GRCh38_STAR_index"
chromsizes="GRCh38_chromSizes"
SalmonIndexDirectory="GRCh38_salmon_index"

wd="multi_cell"

sample_info="processing/rnaseq/multi_cell_line_stm3675/sample_info/multicell_sample_info.csv"
samples=$(cat $sample_info | cut -f1 -d, | tail -n+2 | uniq)

raw_fastq_directory="${wd}/raw_fastq/"
trimmed_fastq_directory="${wd}/trimmed_fastq/"
trimmed_fastqc_directory="${wd}/trimmed_fastqc/"
umi_fastq_directory="${wd}/umi_fastq/"

alignment_directory="${wd}/alignment"
alignment_star_directory="${alignment_directory}/STAR_bam"
alignment_umidedup_bam_directory="${alignment_directory}/umidedup_bam"
alignment_umidedup_bw_directory="${alignment_directory}/umidedup_bw"
alignment_umidedup_fastq_directory="${alignment_directory}/umidedup_fastq"

salmon_counts_directory="${wd}/salmon_counts/"

mkdir -p ${raw_fastq_directory}
mkdir -p ${trimmed_fastq_directory}
mkdir -p ${trimmed_fastqc_directory}
mkdir -p ${umi_fastq_directory}

mkdir -p ${salmon_counts_directory}
mkdir -p ${alignment_directory}
mkdir -p ${alignment_star_directory}
mkdir -p ${alignment_umidedup_bam_directory}
mkdir -p ${alignment_umidedup_bw_directory}
mkdir -p ${alignment_umidedup_fastq_directory}

################################################################################

nextSeqTrim="--nextseq-trim=10"
rd1Adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"

echo $samples
for sample in $samples
do
  echo $sample
  raw_fastq_file="${raw_fastq_directory}/${sample}.fastq.gz"
  ##########################
  echo "processing SE data file: ${sample}.fq.gz"
  ##########################
  echo "extracting UMIs..."
  umi_fastq="${umi_fastq_directory}/${sample}_extractumi.fastq.gz"
  echo $umi_fastq
  zcat $raw_fastq_file | awk 'NR%4==1{ rd_name=$1; rd_info=$2 } NR%4==2{ umi=substr($1,1,10); rd_seq=substr($1,13) } NR%4==0{ print rd_name"_"umi" "rd_info; print rd_seq; print "+"; print substr($1,13) }' | gzip > $umi_fastq
  ##########################
  echo "adapter trimming..."
  trim_fastq="${trimmed_fastq_directory}/${sample}_trim.fastq.gz"
  echo $trim_fastq
  cutadapt -m 20 -O 20 -a "QUALITY=G{20}" -j 8 ${umi_fastq} | cutadapt -m 20 $nextSeqTrim -a $rd1Adapter -j 8 - | cutadapt -m 20 -O 3 -a "r1polyA=A{18}" -j 8 - | cutadapt -m 20 -O 20 -g $rd1Adapter --discard-trimmed -o ${trim_fastq} -j 8 -
  ##########################
  echo "creating quality report..."
  fastqc -o ${trimmed_fastqc_directory} ${trim_fastq}
done

################################################################################

echo $samples
for sample in $samples
do
  echo $sample
  trim_fastq="${trimmed_fastq_directory}/${sample}_trim.fastq.gz"
  ls $trim_fastq
  ##########################
  echo "alignment..."
  align_star="${alignment_star_directory}/${sample}."
  STAR --runThreadN 12 --readFilesCommand zcat --genomeDir $STARIndexDirectory --readFilesIn ${trim_fastq} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.6 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outSAMattributes NH HI NM MD --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ${align_star}
  # umi_tools dedup needs indexed bam files
  samtools index "${align_star}Aligned.sortedByCoord.out.bam"
  ######################################################################
  echo "umi deduplicating..."
  dedup_out="${alignment_umidedup_bam_directory}/${sample}.umidedup.bam"
  dedup_log="${alignment_umidedup_bam_directory}/${sample}.umidedup.log"
  umi_tools dedup -I "${align_star}Aligned.sortedByCoord.out.bam" -S ${dedup_out} --multimapping-detection-method=NH --method=unique --log=${dedup_log}
  samtools index ${dedup_out}
  ######################################################################
  echo "bigwig file generation..."
  bg_pos="${alignment_umidedup_bw_directory}/${sample}.pos.bedGraph"
  bg_neg="${alignment_umidedup_bw_directory}/${sample}.neg.bedGraph"
  bedtools genomecov -ibam $dedup_out -bg -strand + -split > $bg_pos &
  bedtools genomecov -ibam $dedup_out -bg -strand - -split > $bg_neg
  wait
  bg_pos_clip="${alignment_umidedup_bw_directory}/${sample}.pos_clip.bedGraph"
  bg_neg_clip="${alignment_umidedup_bw_directory}/${sample}.neg_clip.bedGraph"
  bedClip $bg_pos $chromsizes $bg_pos_clip &
  bedClip $bg_neg $chromsizes $bg_neg_clip
  wait
  bedSort $bg_pos_clip $bg_pos_clip &
  bedSort $bg_neg_clip $bg_neg_clip
  wait
  bw_pos="${alignment_umidedup_bw_directory}/${sample}.pos.bw"
  bw_neg="${alignment_umidedup_bw_directory}/${sample}.neg.bw"
  bedGraphToBigWig $bg_pos_clip $chromsizes $bw_pos &
  bedGraphToBigWig $bg_neg_clip $chromsizes $bw_neg
  wait
  rm $alignment_umidedup_bw_directory/*.bedGraph
  ######################################################################
  dedup_fastq="${alignment_umidedup_fastq_directory}/${sample}.umidedup.fq.gz"
  #Convert deduplicated bam to randomized fastq for mapping to salmon_index - this takes only primary alignments from bam file (samtools view -F 256) so do not get duplicated reads in fastq when multi-mapped
  samtools fastq $dedup_out | awk '{OFS="\t"; getline seq; getline sep; getline qual; print $0,seq,sep,qual}' | shuf | awk '{OFS="\n"; print $1,$2,$3,$4}' | gzip > $dedup_fastq
done

################################################################################
### Salmon counts

echo $samples
for sample in $samples
do
  echo $sample
  echo "salmon transcript counting..."
  dedup_fq="${alignment_umidedup_fastq_directory}/${sample}.umidedup.fq.gz"
  salmon_counts="${salmon_counts_directory}/${sample}"
  salmon quant -i $SalmonIndexDirectory -l SF -r $dedup_fastq -p 16 -o $salmon_counts
done

################################################################################
###Intron GTF/RO GTF

ir_featurecounts_directory="${wd}/feature_counts/intron_retention/"
ro_featurecounts_directory="${wd}/feature_counts/runon/"
mkdir -p ${ir_featurecounts_directory}
mkdir -p ${ro_featurecounts_directory}

bamfiles=$(ls ${alignment_umidedup_bam_directory}/*bam)
caov3_bamfiles=$(ls ${alignment_umidedup_bam_directory}/*bam | grep -E 'CD|CT')

IntronAnnotationGTF="Hsap_mane_select_introns.saf"

featurecounts_introns="${ir_featurecounts_directory}/featurecounts_introns_multicell"
featureCounts -s 1 -T 8 -F "SAF" -f -O -M -a ${IntronAnnotationGTF} -o ${featurecounts_introns} ${bamfiles}

mane_runon_chunk1000_AnnotationsGTF="Hsap_mane_runon_chunk1000.gtf"
featurecounts_mane_runon_chunk1000="${ro_featurecounts_directory}/featurecounts_mane_runon_chunk1000"

featureCounts -s 1 -T 8 -g "chunk" -O -M -a ${mane_runon_chunk1000_AnnotationsGTF} -o ${featurecounts_mane_runon_chunk1000} ${bamfiles}

mane_runon_joined_chunk_above_caov3_nc10_AnnotationsGTF="Hsap_mane_runon_joined_chunk1000_above_caov3_nc10.gtf"
featurecounts_mane_runon_joined_chunk_above_caov3_nc10="${ro_featurecounts_directory}/featurecounts_mane_runon_joined_chunk_above_caov3_nc10"

featureCounts -s 1 -T 8 -O -M -a ${mane_runon_joined_chunk_above_caov3_nc10_AnnotationsGTF} -o ${featurecounts_mane_runon_joined_chunk_above_caov3_nc10} ${caov3_bamfiles}

################################################################################

mane_runon_joined_chunk_above_allcells_nc10_AnnotationsGTF="Hsap_mane_runon_joined_chunk1000_above_allcells_nc10.gtf"
featurecounts_mane_runon_joined_chunk_above_allcells_nc10="${ro_featurecounts_directory}/featurecounts_mane_runon_joined_chunk_above_allcells_nc10"

featureCounts -s 1 -T 8 -O -M -a ${mane_runon_joined_chunk_above_allcells_nc10_AnnotationsGTF} -o ${featurecounts_mane_runon_joined_chunk_above_allcells_nc10} ${bamfiles}
