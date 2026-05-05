#!/usr/bin/env bash
set -euo pipefail

sample="$1"
trimmed_fastq="$2"
out_dir="$3"

mkdir -p "${out_dir}"
atog_fastq="${out_dir}/${sample}_AtoG.fq"
backcsv="${out_dir}/${sample}_GtoA_backconversion.csv"
basecount="${out_dir}/${sample}_trim_basecount.csv"

rm -f "${atog_fastq}" "${backcsv}" "${basecount}"

zcat "${trimmed_fastq}" | awk -v basecount="${basecount}" -v atogfq="${atog_fastq}" -v gtoacsv="${backcsv}" 'BEGIN{print "A,C,G,T,N" > basecount} NR%4==1{ rd_name=$0; mod_rd_name=gensub("@", "", 1, $1) } NR%4==2{ seq=$1; count_seq=$1; mod_seq=$1; acount=gsub("A","",count_seq); ccount=gsub("C","",count_seq); gcount=gsub("G","",count_seq); tcount=gsub("T","",count_seq); ncount=gsub("N","",count_seq); print acount","ccount","gcount","tcount","ncount >> basecount } NR%4==0{ if (gsub("A","G",mod_seq) < 4) {print rd_name >> atogfq; print mod_seq >> atogfq; print "+" >> atogfq; print $0 >> atogfq; print mod_rd_name "," seq >> gtoacsv}}'

gzip -c "${atog_fastq}" > "${atog_fastq}.gz"
rm -f "${atog_fastq}"
