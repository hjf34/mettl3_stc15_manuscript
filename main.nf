#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.sample_info = params.sample_info ?: 'processing/glori/sample_info/glori_sample_info.csv'
params.combined_fastq_dir = params.combined_fastq_dir ?: 'processing/glori/combined_fastq'
params.trimmed_fastq_dir = params.trimmed_fastq_dir ?: 'processing/glori/trimmed_fastq'
params.atog_fastq_dir = params.atog_fastq_dir ?: 'processing/glori/atog_fastq'
params.alignment_dir = params.alignment_dir ?: 'processing/glori/alignment'
params.si_alignment_dir = params.si_alignment_dir ?: "${params.alignment_dir}/spikein"
params.genome_alignment_dir = params.genome_alignment_dir ?: "${params.alignment_dir}/genome"
params.si_alignment_processing_dir = params.si_alignment_processing_dir ?: "${params.si_alignment_dir}/processing"
params.genome_alignment_processing_dir = params.genome_alignment_processing_dir ?: "${params.genome_alignment_dir}/processing"

params.genome_fasta = params.genome_fasta ?: 'processing/glori/genome_setup/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz'
params.genome_gtf = params.genome_gtf ?: 'processing/glori/genome_setup/Homo_sapiens.GRCh38.108.gtf'
params.combined_genome_fasta = params.combined_genome_fasta ?: 'processing/glori/genome_setup/GRCh38_DNA_combined_strands_AtoG.fa'
params.combined_annotations = params.combined_annotations ?: 'processing/glori/genome_setup/GRCh38_annotations_combined_positive_and_negative_revcomp.gtf'
params.genome_star_index = params.genome_star_index ?: 'processing/glori/genome_setup/STARIndex_combined_strands_AtoG_withannot'
params.si_atog_fasta = params.si_atog_fasta ?: 'processing/glori/resources/spike_in/atog_spikein.fa'
params.si_star_index = params.si_star_index ?: 'processing/glori/genome_setup/STARIndex_spike_in'
params.mpileup_dir = params.mpileup_dir ?: 'processing/glori/mpileup'
params.genome_mpileup_dir = params.genome_mpileup_dir ?: "${params.mpileup_dir}/genome"
params.si_mpileup_dir = params.si_mpileup_dir ?: "${params.mpileup_dir}/spikein"
params.chunks_dir = params.chunks_dir ?: 'processing/glori/genome_setup/chunks100'
params.root = params.root ?: projectDir

workflow {
    samples = Channel
        .fromPath(params.sample_info)
        .flatMap { file ->
            file.readLines().drop(1).collect { line -> line.split(',') }
        }
        .map { tokens -> tuple(tokens[0], tokens) }

    indexReady = generateGenomeIndex(params.genome_fasta, params.genome_gtf, params.combined_genome_fasta, params.combined_annotations, params.genome_star_index)
        .combine(generateSpikeInIndex(params.si_atog_fasta, params.si_star_index))
        .map { 'ready' }
        .broadcast()

    trimmed = trimFastq(samples, params.combined_fastq_dir, params.trimmed_fastq_dir)
    atog = convertFastqAtoG(trimmed, params.atog_fastq_dir)

    alignSpikeIn(indexReady, atog, params.si_star_index, params.si_alignment_dir)
    genomeAlign = alignGenome(indexReady, atog, params.genome_star_index, params.genome_alignment_dir)
    backConvertSpikeIn(alignSpikeIn)
    backConvertGenome(genomeAlign)
    spikeInMpileup(backConvertSpikeIn)
    runMpileup(backConvertGenome.out.collect())
    combineMpileupChunks(runMpileup.out)
}

process generateGenomeIndex {
    tag 'glori-genome-index'
    cpus 8
    memory '32 GB'

    input:
    path genome_fasta
    path genome_gtf

    output:
    val('ready')

    script:
    def root = params.root
    def combined_fasta = "${root}/${params.combined_genome_fasta}"
    def combined_gtf = "${root}/${params.combined_annotations}"
    def index_dir = "${root}/${params.genome_star_index}"
    def chunks_dir = "${root}/${params.chunks_dir}"
    """
    mkdir -p ${index_dir}
    mkdir -p ${chunks_dir}
    Rscript ${root}/processing/glori/genome_setup/01_GRCh38_genome_reverse_complement_and_AtoG_conversion.R
    Rscript ${root}/processing/glori/genome_setup/03_GRCh38_genome_positive_strand_and_negative_reverse_complement_strand_transcriptome_A_positions.R
    bash ${root}/processing/glori/genome_setup/04_split_A_positions_into_100_line_chunks_for_parallel_mpileup.sh
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} \
         --genomeDir ${index_dir} \
         --genomeFastaFiles ${combined_fasta} \
         --sjdbGTFfile ${combined_gtf} \
         --sjdbOverhang 100 \
         --genomeSuffixLengthMax 300
    """
}

process generateSpikeInIndex {
    tag 'glori-spikein-index'
    cpus 4
    memory '16 GB'

    input:
    path si_atog_fasta

    output:
    val('ready')

    script:
    def root = params.root
    def index_dir = "${root}/${params.si_star_index}"
    """
    mkdir -p ${index_dir}
    STAR --runMode genomeGenerate --runThreadN ${task.cpus} \
         --genomeDir ${index_dir} \
         --genomeFastaFiles ${si_atog_fasta} \
         --genomeSAindexNbases 3
    """
}

process trimFastq {
    tag { sample }
    cpus 8
    memory '16 GB'

    input:
    tuple val(sample), val(meta)

    output:
    tuple val(sample), path("${params.root}/${params.trimmed_fastq_dir}/${sample}_trim.fq.gz"), path("${params.root}/${params.trimmed_fastq_dir}/${sample}_trim_basecount.csv")

    script:
    def root = params.root
    def input_fastq = "${root}/${params.combined_fastq_dir}/${sample}.fq.gz"
    def output_dir = "${root}/${params.trimmed_fastq_dir}"
    """
    mkdir -p ${output_dir}
    cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 10 -m 20 -j ${task.cpus} --trimmed-only \
        -o ${output_dir}/${sample}_trim.fq.gz ${input_fastq}
    """
}

process convertFastqAtoG {
    tag { sample }
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample), path(trimmed_fastq), path(trimmed_count)

    output:
    tuple val(sample), path("${params.root}/${params.atog_fastq_dir}/${sample}_AtoG.fq.gz"), path("${params.root}/${params.atog_fastq_dir}/${sample}_GtoA_backconversion.csv")

    script:
    def root = params.root
    def output_dir = "${root}/${params.atog_fastq_dir}"
    """
    mkdir -p ${output_dir}
    ${root}/bin/convert_fastq_to_atog.sh ${sample} ${trimmed_fastq} ${output_dir}
    """
}

process alignSpikeIn {
    tag { sample }
    cpus 16
    memory '32 GB'

    input:
    val(dummy) from indexReady
    tuple val(sample), path(atog_fastq), path(backconversion_csv)

    output:
    tuple val(sample), path("${params.root}/${params.si_alignment_dir}/${sample}_AtoG_alignSI.Aligned.sortedByCoord.out.bam"), path(backconversion_csv)

    script:
    def root = params.root
    def index_dir = "${root}/${params.si_star_index}"
    def out_dir = "${root}/${params.si_alignment_dir}"
    """
    mkdir -p ${out_dir}
    STAR --runThreadN ${task.cpus} --readFilesCommand zcat \
         --genomeDir ${index_dir} \
         --readFilesIn ${atog_fastq} \
         --outFilterType BySJout \
         --outFilterMultimapNmax 1 \
         --alignSJoverhangMin 999 \
         --alignSJDBoverhangMin 999 \
         --outFilterMismatchNmax 2 \
         --outFilterMismatchNoverLmax 0.3 \
         --alignIntronMin 1 \
         --alignIntronMax 1 \
         --outSAMtype BAM SortedByCoordinate \
         --alignEndsType EndToEnd \
         --limitBAMsortRAM 1027490585 \
         --outFileNamePrefix ${out_dir}/${sample}_AtoG_alignSI.
    """
}

process alignGenome {
    tag { sample }
    cpus 16
    memory '32 GB'

    input:
    val(dummy) from indexReady
    tuple val(sample), path(atog_fastq), path(backconversion_csv)

    output:
    tuple val(sample), path("${params.root}/${params.genome_alignment_dir}/${sample}_AtoG_alignGenome.Aligned.out.bam"), path(backconversion_csv)

    script:
    def root = params.root
    def index_dir = "${root}/${params.genome_star_index}"
    def out_dir = "${root}/${params.genome_alignment_dir}"
    """
    mkdir -p ${out_dir}
    STAR --runThreadN ${task.cpus} --readFilesCommand zcat \
         --genomeDir ${index_dir} \
         --readFilesIn ${atog_fastq} \
         --outFilterType BySJout \
         --outFilterMultimapNmax 10 \
         --outFilterMismatchNmax 2 \
         --alignIntronMin 20 \
         --alignIntronMax 1000000 \
         --outMultimapperOrder Random \
         --outSAMattributes All \
         --outSAMtype BAM Unsorted \
         --outFileNamePrefix ${out_dir}/${sample}_AtoG_alignGenome.
    """
}

process backConvertGenome {
    tag { sample }
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample), path(genome_bam), path(backconversion_csv)

    output:
    path("${params.root}/${params.genome_alignment_processing_dir}/${sample}_alignGenome_gtoa.bam"), path("${params.root}/${params.genome_alignment_processing_dir}/${sample}_alignGenome_gtoa.bam.bai")

    script:
    def root = params.root
    def output_dir = "${root}/${params.genome_alignment_processing_dir}"
    """
    mkdir -p ${output_dir}
    ${root}/bin/back_convert_bam.sh ${sample} genome ${genome_bam} ${backconversion_csv} ${output_dir}
    """
}

process backConvertSpikeIn {
    tag { sample }
    cpus 4
    memory '16 GB'

    input:
    tuple val(sample), path(si_bam), path(backconversion_csv)

    output:
    tuple val(sample), path("${params.root}/${params.si_alignment_processing_dir}/${sample}_alignSI_gtoa.bam"), path("${params.root}/${params.si_alignment_processing_dir}/${sample}_alignSI_gtoa.bam.bai")

    script:
    def root = params.root
    def output_dir = "${root}/${params.si_alignment_processing_dir}"
    """
    mkdir -p ${output_dir}
    ${root}/bin/back_convert_bam.sh ${sample} spikein ${si_bam} ${backconversion_csv} ${output_dir}
    """
}

process spikeInMpileup {
    tag { sample }
    cpus 4
    memory '8 GB'

    input:
    tuple val(sample), path(si_bam), path(si_bai)

    output:
    path("${params.root}/${params.si_mpileup_dir}/${sample}_alignSI_all_mpileup_allpositions.txt")

    script:
    def root = params.root
    def output_file = "${root}/${params.si_mpileup_dir}/${sample}_alignSI_all_mpileup_allpositions.txt"
    def si_fasta = "${root}/${params.si_atog_fasta}"
    """
    mkdir -p ${root}/${params.si_mpileup_dir}
    bcftools mpileup -f ${si_fasta} -d 1000000 -O u --annotate FORMAT/DP,FORMAT/AD ${si_bam} | bcftools query -i 'FORMAT/DP>0' -f '%CHROM\t%POS\t%REF\t%ALT[\t%DP\t%AD]\n' > ${output_file}
    """
}

process runMpileup {
    tag 'glori-mpileup'
    cpus 40
    memory '64 GB'

    input:
    path(bam_files)

    output:
    path("${params.root}/${params.genome_mpileup_dir}/allsamples_alignGenome_chunks100_mpileup/*.tsv")

    script:
    def root = params.root
    def bam_list = "all_alignGenome_gtoa_bamfiles.txt"
    def out_dir = "${root}/${params.genome_mpileup_dir}/allsamples_alignGenome_chunks100_mpileup"
    def ref = "${root}/${params.combined_genome_fasta}"
    """
    mkdir -p ${out_dir}
    ls ${root}/${params.genome_alignment_processing_dir}/*_alignGenome_gtoa.bam > ${bam_list}
    cd ${root}/processing/glori
    REF="${ref}" BAM="${bam_list}" CHUNK_DIR="genome_setup/chunks100" OUTDIR="${out_dir}" bash mpileup_at_transcriptome_A_positions_in_100_line_chunks.sh ${task.cpus}
    """
}

process combineMpileupChunks {
    tag 'combine-mpileup-chunks'
    cpus 4
    memory '16 GB'

    input:
    path(mpileup_chunks)

    output:
    path("${params.root}/${params.genome_mpileup_dir}/allsamples_alignGenome_chunks100_mpileup/combined_*_chunks.tsv")

    script:
    def root = params.root
    def out_dir = "${root}/${params.genome_mpileup_dir}/allsamples_alignGenome_chunks100_mpileup"
    """
    cd ${out_dir}
    OUTDIR="${out_dir}" bash ${root}/processing/glori/02_combine_all_mpileup_100_line_chunks_per_chromosome.sh
    """
}
