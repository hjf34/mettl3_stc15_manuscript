################################################################################
#### Generate STAR index for Human genome GRCh38

# Download GRCh38 genome and annotation file from Ensembl https://ftp.ensembl.org/
STARIndexDirectory="GRCh38_STAR_index"
GenomeFasta="Homo_sapiens.GRCh38.dna.primary_assembly.fa"
GenomeFai="Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai"
AnnotationsGTF="Homo_sapiens.GRCh38.108.gtf"
TranscriptomeFasta="GRCh38_transcriptome.fa"

#Index genome fasta
samtools faidx $GenomeFasta
#Generate chromosome sizes file
cut -f1,2 $GenomeFai > "GRCh38_chromSizes"

#Generate STAR alignment index
STAR --runMode genomeGenerate --runThreadN 8 --genomeDir $STARIndexDirectory --genomeFastaFiles $GenomeFasta --sjdbGTFfile $AnnotationsGTF --sjdbOverhang 100

################################################################################
#### Generate Salmon index for Human genome GRCh38

gffread -w $TranscriptomeFasta -g $GenomeFasta $AnnotationsGTF
awk -i inplace '/^>/{print $1; next}{print}' $TranscriptomeFasta

DecoyFile="decoys.txt"
GentromeFile="gentrome.fa"
SalmonIndexDirectory="GRCh38_salmon_index"
grep "^>" $GenomeFasta | cut -d " " -f 1 > $DecoyFile
sed -i.bak -e 's/>//g' $DecoyFile
cat $TranscriptomeFasta $GenomeFasta > $GentromeFile
salmon index -t $GentromeFile -d $DecoyFile -p 12 -i $SalmonIndexDirectory
