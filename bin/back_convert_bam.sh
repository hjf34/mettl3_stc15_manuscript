#!/usr/bin/env bash
set -euo pipefail

sample="$1"
tag="$2"
bam="$3"
gtoa_csv="$4"
out_dir="$5"

mkdir -p "${out_dir}"
tmp_sam="${out_dir}/${sample}_${tag}_positive_only.sam"
tmp_names="${out_dir}/${sample}_${tag}_positive_only_names.txt"
tmp_all="${out_dir}/${sample}_${tag}_allsequences_gtoa.txt"
tmp_gtoa_sam="${out_dir}/${sample}_${tag}_gtoa.sam"
out_bam="${out_dir}/${sample}_${tag}_gtoa.bam"

samtools view -F 16 "${bam}" > "${tmp_sam}"
cut -f1 "${tmp_sam}" > "${tmp_names}"
awk -F'[,]' 'NR==FNR{a[$1]=$2; next} {print $1, a[$1]}' "${gtoa_csv}" "${tmp_names}" > "${tmp_all}"
samtools view -H "${bam}" > "${tmp_gtoa_sam}"
awk 'FNR==NR{a[NR]=$2;next}{$10=a[FNR]}1' OFS='\t' "${tmp_all}" "${tmp_sam}" >> "${tmp_gtoa_sam}"
samtools sort "${tmp_gtoa_sam}" -o "${out_bam}"
samtools index "${out_bam}"
rm -f "${tmp_sam}" "${tmp_names}" "${tmp_all}" "${tmp_gtoa_sam}"
