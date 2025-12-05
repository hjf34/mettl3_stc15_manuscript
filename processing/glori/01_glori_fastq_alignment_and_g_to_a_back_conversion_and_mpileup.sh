################################################################################
#### Directory set up

wd=""

sample_info="sample_info.csv"
samples=$(cat $sample_info | cut -f1 -d, | tail -n+2 | uniq)

combined_fastq_directory="${wd}/combined_fastq/"
trimmed_fastq_directory="${wd}/trimmed_fastq/"
trimmed_fastqc_directory="${wd}/trimmed_fastqc/"
atog_fastq_directory="${wd}/atog_fastq/"

alignment_directory="${wd}/alignment/"
si_alignment_directory="${alignment_directory}/spikein/"
genome_alignment_directory="${alignment_directory}/genome/"
si_alignment_processing_directory="${alignment_directory}/spikein/processing/"
genome_alignment_processing_directory="${alignment_directory}/genome/processing/"

mpileup_directory="${wd}/mpileup/"
si_mpileup_directory="${mpileup_directory}/spikein/"
genome_mpileup_directory="${mpileup_directory}/genome/"

mkdir -p ${combined_fastq_directory}
mkdir -p ${trimmed_fastq_directory}
mkdir -p ${trimmed_fastqc_directory}
mkdir -p ${atog_fastq_directory}

mkdir -p ${alignment_directory}
mkdir -p ${si_alignment_directory}
mkdir -p ${genome_alignment_directory}
mkdir -p ${si_alignment_processing_directory}
mkdir -p ${genome_alignment_processing_directory}

mkdir -p ${mpileup_directory}
mkdir -p ${si_mpileup_directory}
mkdir -p ${genome_mpileup_directory}

# Raw fastq files in combined_fastq_directory (concatenate fastq files from multiple lanes). Use read1 of PE150 sequencing reads

################################################################################
#### Trimming adapters, minimum overlap of 10 bases between read and adapter, discard reads shorter than 20 bases after trimming or untrimmed reads

r1adapter="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC"
echo $samples
for sample in $samples
do
  echo $sample
  combined_fastq="${combined_fastq_directory}/${sample}.fq.gz"
  echo $combined_fastq
  trimmed_fastq="${trimmed_fastq_directory}/${sample}_trim.fq.gz"
  echo $trimmed_fastq
  trimmed_report="${trimmed_fastq_directory}/${sample}_trim.json"
  cutadapt -a ${r1adapter} -O 10 -m 20 -j 16 --trimmed-only -o ${trimmed_fastq} --json=${trimmed_report} ${combined_fastq}
done

################################################################################
#### FastQC on trimmed files
trimmed_fastq_files=$(ls ${trimmed_fastq_directory}/*.gz)
echo $trimmed_fastq_files
fastqc -o ${trimmed_fastqc_directory} ${trimmed_fastq_files}

################################################################################
#### Generate stats on ACGTN read count in each trimmed fastq file
#### Filter reads with more than 3 As
#### Convert all As to Gs in reads in fastq file
#### Generate csv file for converting Gs back to As in aligned reads

echo $samples
for sample in $samples
do
  ##############################################################################
  echo $sample
  trimmed_fastq="${trimmed_fastq_directory}/${sample}_trim.fq.gz"
  echo $trimmed_fastq
  ##############################################################################
  trimmed_fastq_basecount="${trimmed_fastq_directory}/${sample}_trim_basecount.csv"
  echo $trimmed_fastq_basecount
  atog_fastq="${atog_fastq_directory}/${sample}_AtoG.fq"
  echo $atog_fastq
  gtoa_backconversion_csv="${atog_fastq_directory}/${sample}_GtoA_backconversion.csv"
  echo $gtoa_backconversion_csv
  rm $trimmed_fastq_basecount
  rm $atog_fastq
  rm $gtoa_backconversion_csv
  ##############################################################################
  zcat $trimmed_fastq | awk -v basecount=$trimmed_fastq_basecount -v atogfq=$atog_fastq -v gtoacsv=$gtoa_backconversion_csv 'BEGIN{print "A,C,G,T,N" > basecount} NR%4==1{ rd_name=$0; mod_rd_name=gensub("@", "", 1, $1) } NR%4==2{ seq=$1; count_seq=$1; mod_seq=$1; acount=gsub("A","",count_seq); ccount=gsub("C","",count_seq); gcount=gsub("G","",count_seq); tcount=gsub("T","",count_seq); ncount=gsub("N","",count_seq); print acount","ccount","gcount","tcount","ncount >> basecount} NR%4==0{ if (gsub("A","G",mod_seq) < 4) {print rd_name >> atogfq; print mod_seq >> atogfq; print "+" >> atogfq; print $0 >> atogfq; print mod_rd_name","seq >> gtoacsv}}'
  gzip $atog_fastq
done

################################################################################
#### Alignment to A to G converted spike-in sequences
#### Alignment to A to G converted GRCh38 genome (A to G convert both the positive and reverse complemented negative strand then combine strands into one genome STAR index with converted gtf annotation file)

si_star_index="STARIndex_spike_in"
genome_star_index="STARIndex_combined_strands_AtoG_withannot"

samples=$(cat $sample_info | cut -f1 -d, | uniq)
echo $samples
for sample in $samples
do
  ##############################################################################
  echo $sample
  atog_fastq="${atog_fastq_directory}/${sample}_AtoG.fq.gz"
  echo $atog_fastq
  align_si="${si_alignment_directory}/${sample}_AtoG_alignSI."
  align_genome="${genome_alignment_directory}/${sample}_AtoG_alignGenome."
  ##############################################################################
  ##Align to spike-in sequences
  STAR --runThreadN 16 --readFilesCommand zcat --genomeDir ${si_star_index} --readFilesIn ${atog_fastq} --outFilterType BySJout --outFilterMultimapNmax 1 --alignSJoverhangMin 999 --alignSJDBoverhangMin 999 --outFilterMismatchNmax 2 --outFilterMismatchNoverLmax 0.3 --alignIntronMin 1 --alignIntronMax 1 --outSAMtype BAM SortedByCoordinate --alignEndsType EndToEnd --limitBAMsortRAM 1027490585 --outFileNamePrefix ${align_si}
  ##############################################################################
  ##Align to genome
  STAR --runThreadN 16 --readFilesCommand zcat --genomeDir ${genome_star_index} --readFilesIn ${atog_fastq} --outFilterType BySJout --outFilterMultimapNmax 10 --outFilterMismatchNmax 2 --alignIntronMin 20 --alignIntronMax 1000000 --outMultimapperOrder Random --outSAMattributes All --outSAMtype BAM Unsorted --outFileNamePrefix ${align_genome}
done

################################################################################
#### Filter out reads aligned to negative strand (most reads should align in a positive orientation either to the positive strand or the reverse complemented negative strand)
#### Back convert Gs back to As

echo $samples
for sample in $samples
do
  ##############################################################################
  echo $sample
  align_si_bam="${si_alignment_directory}/${sample}_AtoG_alignSI.Aligned.sortedByCoord.out.bam"
  align_genome_bam="${genome_alignment_directory}/${sample}_AtoG_alignGenome.Aligned.out.bam"
  gtoa_backconversion_csv="${atog_fastq_directory}/${sample}_GtoA_backconversion.csv"
  ls $align_si_bam
  ls $align_genome_bam
  ls $gtoa_backconversion_csv
  ##############################################################################
  align_si_positive_only_sam="${si_alignment_processing_directory}/${sample}_alignSI_positive_only.sam"
  align_si_positive_only_names="${si_alignment_processing_directory}/${sample}_alignSI_positive_only_names.txt"
  align_si_allsequences_gtoa="${si_alignment_processing_directory}/${sample}_alignSI_allsequences_gtoa.txt"
  align_si_gtoa_sam="${si_alignment_processing_directory}/${sample}_alignSI_gtoa.sam"
  align_si_gtoa_bam="${si_alignment_processing_directory}/${sample}_alignSI_gtoa.bam"
  ##############################################################################
  align_genome_positive_only_sam="${genome_alignment_processing_directory}/${sample}_alignGenome_positive_only.sam"
  align_genome_positive_only_names="${genome_alignment_processing_directory}/${sample}_alignGenome_positive_only_names.txt"
  align_genome_allsequences_gtoa="${genome_alignment_processing_directory}/${sample}_alignGenome_allsequences_gtoa.txt"
  align_genome_gtoa_sam="${genome_alignment_processing_directory}/${sample}_alignGenome_gtoa.sam"
  align_genome_gtoa_bam="${genome_alignment_processing_directory}/${sample}_alignGenome_gtoa.bam"
  ##############################################################################
  samtools view -F 16 ${align_si_bam} > ${align_si_positive_only_sam}
  cut -d$'\t' -f1 ${align_si_positive_only_sam} > ${align_si_positive_only_names}
  ##############################################################################
  samtools view -F 16 ${align_genome_bam} > ${align_genome_positive_only_sam}
  cut -d$'\t' -f1 ${align_genome_positive_only_sam} > ${align_genome_positive_only_names}
  ##############################################################################
  awk -F'[,]' 'NR==FNR{a[$1]=$2; next} {print $1, a[$1]}' ${gtoa_backconversion_csv} ${align_si_positive_only_names} > ${align_si_allsequences_gtoa}
  awk -F'[,]' 'NR==FNR{a[$1]=$2; next} {print $1, a[$1]}' ${gtoa_backconversion_csv} ${align_genome_positive_only_names} > ${align_genome_allsequences_gtoa}
  ##############################################################################
  samtools view -H ${align_si_bam} > ${align_si_gtoa_sam}
  awk 'FNR==NR{a[NR]=$2;next}{$10=a[FNR]}1' OFS='\t' ${align_si_allsequences_gtoa} ${align_si_positive_only_sam} >> ${align_si_gtoa_sam}
  samtools sort ${align_si_gtoa_sam} -o ${align_si_gtoa_bam}
  samtools index ${align_si_gtoa_bam}
  ##############################################################################
  samtools view -H ${align_genome_bam} > ${align_genome_gtoa_sam}
  awk 'FNR==NR{a[NR]=$2;next}{$10=a[FNR]}1' OFS='\t' ${align_genome_allsequences_gtoa} ${align_genome_positive_only_sam} >> ${align_genome_gtoa_sam}
  samtools sort ${align_genome_gtoa_sam} -o ${align_genome_gtoa_bam}
  samtools index ${align_genome_gtoa_bam}
  ##############################################################################
  rm $align_si_positive_only_sam $align_si_positive_only_names $align_si_allsequences_gtoa $align_si_gtoa_sam
  rm $align_genome_positive_only_sam $align_genome_positive_only_names $align_genome_allsequences_gtoa $align_genome_gtoa_sam
done

ls -d ${genome_alignment_processing_directory}/*bam > "${genome_alignment_processing_directory}/all_alignGenome_gtoa_bamfiles.txt"

################################################################################
#### Counting As (m6As) and Gs (unmodified As) at A positions in spike in and genome aligned bam files

################################################################################
#### bcftools mpileup at all spike in sequences
################################################################################
###Call number of As and Gs at each position using mpileup

spikein_atog_fasta="processing/glori/resources/spike_in/atog_spikein.fa"

###Spike-in mpileup

samples=$(cat $sample_info | cut -f1 -d, | uniq)
echo $samples
for sample in $samples
do
  ##############################################################################
  echo $sample
  align_si_gtoa_bam="${si_alignment_processing_directory}/${sample}_alignSI_gtoa.bam"
  ls ${align_si_gtoa_bam}
  mpileup_si_all="${si_mpileup_directory}/${sample}_alignSI_all_mpileup_allpositions.txt"
  echo ${mpileup_si_all}
  bcftools mpileup -f ${spikein_atog_fasta} -d 1000000 -O u --annotate FORMAT/DP,FORMAT/AD ${align_si_gtoa_bam} | bcftools query -i 'FORMAT/DP>0' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' > ${mpileup_si_all}
done

################################################################################
#### bcftools mpileup at transcriptome A positions in chunks of 100 to allow efficient parallelisation

bash mpileup_at_transcriptome_A_positions_in_100_line_chunks.sh 40
