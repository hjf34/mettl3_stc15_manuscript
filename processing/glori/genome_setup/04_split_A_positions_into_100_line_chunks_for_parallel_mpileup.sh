# Process both chrn and chrp files
for file in apos_exons_chrn_*.txt apos_exons_chrp_*.txt; do
    [ -e "$file" ] || continue

    echo "Processing $file..."

    base="${file%.txt}"
    chunk_num=0
    line_count=0
    chunk_lines=()

    while IFS= read -r line; do
        chunk_lines+=("$line")
        ((line_count++))

        if (( line_count == 100 )); then
            printf "%s\n" "${chunk_lines[@]}" > "${base}_chunk_$(printf "%06d" $chunk_num).txt"
            chunk_lines=()
            line_count=0
            ((chunk_num++))
        fi
    done < "$file"

    # Write final chunk if it has leftover lines
    if (( line_count > 0 )); then
        printf "%s\n" "${chunk_lines[@]}" > "${base}_chunk_$(printf "%06d" $chunk_num).txt"
    fi

    echo "Finished $file → $(($chunk_num + 1)) chunks"
done
