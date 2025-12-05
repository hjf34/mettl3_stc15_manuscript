chroms=$(ls apos_exons_chr*_*_chunk_*.tsv | sed -E 's/^apos_exons_(chr[^_]+_[^_]+)_chunk_.*$/_\1_/' | sort -u)

# Define chromosome prefixes to group by
for prefix in $chroms; do
  echo $prefix

  OUTPUT_FILE="$OUTDIR/combined${prefix}chunks.tsv"
  > "$OUTPUT_FILE"

  # Find and combine matching files
  find "$OUTDIR" -maxdepth 1 -name "apos_exons${prefix}chunk*.tsv" | sort | while read -r file; do
      cat "$file" >> "$OUTPUT_FILE"
  done

  echo "✅ Combined $prefix chunks into $OUTPUT_FILE"
done
