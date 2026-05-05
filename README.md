# METTL3 STC-15 Manuscript Scripts

## Storm Therapeutics

### Harry Fischl

---

Scripts to process and analyse GLORI and RNA-seq data for the manuscript:
"Inhibition of METTL3 by STC-15 induces RNA misprocessing that results in dsRNA formation and activates innate immunity"

This repository also contains a Nextflow pipeline to process GLORI data

---

### GLORI pipeline details

- GLORI processing steps:
  - generate STAR index from A-to-G-converted, combined positive and reverse-complemented genome strands
  - generate STAR index from A-to-G-converted spike in probes
  - fastq trimming (cutadapt)
  - fastq A to G conversion and removal of reads with >3 As
  - fastq alignment to A to G converted genome (STAR)
  - fastq alignment to A to G converted spike-in sequences (STAR)
  - back conversion of aligned reads G to A
  - count As and Gs at each A position in transcriptome (parallelized bcftools mpileup)
  - calculate m6A proportions A/(A+G)
  - filter for positions with A+G count ≥ 20 and m6A proportion ≥ 0.25

### Software Versions

- STAR v2.7.10a
- cutadapt v4.8
- bcftools v1.17
- samtools v1.17

---

### Nextflow pipeline
A Nextflow workflow has been added for GLORI genome generation, fastq A-to-G conversion, STAR spike-in alignment, genome alignment, and back-conversion. Run from the repository root with:

```bash
nextflow run main.nf
```

---

Raw data: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312944