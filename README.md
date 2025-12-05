# METTL3 STC-15 Manuscript Scripts

Scripts to process and analyse GLORI and RNA-seq data for the manuscript:
"Inhibition of METTL3 by STC-15 induces RNA misprocessing that results in dsRNA formation and activates innate immunity"

---

## Storm Therapeutics

### Harry Fischl

---

- GLORI processing steps:
  - generate STAR index from A-to-G-converted, combined positive and reverse-complemented genome strands
  - fastq A to G conversion
  - fastq alignment to A to G converted genome (STAR)
  - back conversion of aligned reads G to A
  - count As and Gs at each A position in transcriptome (parallelized bcftools mpileup)
  - calculate m6A proportions A/A+G
  - filter for positions with A+G count ≥ 20 and m6A proportion ≥ 0.25
