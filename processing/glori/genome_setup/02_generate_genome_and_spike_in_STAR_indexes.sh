################################################################################
#### STAR genome generate on A to G converted genome
#### A to G convert both the positive and reverse complemented negative strand then combine strands into one genome STAR index with converted gtf annotation file

genome_star_index="STARIndex_combined_strands_AtoG_withannot"
combined_strands_genome_AtoG_fasta="GRCh38_DNA_combined_strands_AtoG.fa"
combined_annotations="GRCh38_annotations_combined_positive_and_negative_revcomp.gtf"

#### -takes ~3 days
STAR --runMode genomeGenerate --runThreadN 32 --genomeDir $genome_star_index --genomeFastaFiles $combined_strands_genome_AtoG_fasta --sjdbGTFfile $combined_annotations --sjdbOverhang 100 --genomeSuffixLengthMax 300

#### STAR genome generate on A to G converted spike in sequences
si_star_index="STARIndex_spike_in"
spikein_atog_fasta="processing/glori/resources/spike_in/atog_spikein.fa"

STAR --runMode genomeGenerate --runThreadN 16 --genomeDir $si_star_index --genomeFastaFiles $spikein_atog_fasta --genomeSAindexNbases 3
