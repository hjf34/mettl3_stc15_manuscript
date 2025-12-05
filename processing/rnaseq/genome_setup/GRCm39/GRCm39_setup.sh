################################################################################
#### Generate STAR index for Mouse genome GRCm39

# Download GRCm39 genome and annotation file from Ensembl https://ftp.ensembl.org/
STARIndexDirectory="GRCm39_STAR_index"
GenomeFasta="Mus_musculus.GRCm39.dna.primary_assembly.fa"
GenomeFai="Mus_musculus.GRCm39.dna.primary_assembly.fa.fai"
AnnotationsGTF="Mus_musculus.GRCm39.111.gtf"
TranscriptomeFasta="GRCm39_transcriptome.fa"

#Index genome fasta
samtools faidx $GenomeFasta
#Generate chromosome sizes file
cut -f1,2 $GenomeFai > "GRCm39_chromSizes"

#Generate STAR alignment index
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $STARIndexDirectory --genomeFastaFiles $GenomeFasta --sjdbGTFfile $AnnotationsGTF --sjdbOverhang 100

################################################################################
#### Generate Salmon index for Mouse genome GRCm39

gffread -w $TranscriptomeFasta -g $GenomeFasta $AnnotationsGTF
awk -i inplace '/^>/{print $1; next}{print}' $TranscriptomeFasta

DecoyFile="decoys.txt"
GentromeFile="gentrome.fa"
SalmonIndexDirectory="GRCm39_salmon_index"
grep "^>" $GenomeFasta | cut -d " " -f 1 > $DecoyFile
sed -i.bak -e 's/>//g' $DecoyFile
cat $TranscriptomeFasta $GenomeFasta > $GentromeFile
salmon index -t $GentromeFile -d $DecoyFile -p 12 -i $SalmonIndexDirectory
