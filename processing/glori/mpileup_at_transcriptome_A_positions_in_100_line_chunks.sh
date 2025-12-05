#!/bin/bash

# === USAGE CHECK ===
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <number_of_parallel_jobs>"
    exit 1
fi

JOBS="$1"

# === CONFIGURATION ===
REF="GRCh38_DNA_combined_strands_AtoG.fa"
BAM="all_alignGenome_gtoa_bamfiles.txt"
CHUNK_DIR="chunks100" # A positions in 100 line chunks by chromosome
OUTDIR="mpileup/genome/allsamples_alignGenome_chunks100_mpileup"
mkdir -p $OUTDIR

# === FUNCTION TO RUN MPILEUP ===
run_mpileup() {
    region_file="$1"
    chunk_name=$(basename "$region_file" .txt)
    output_file="$OUTDIR/${chunk_name}.tsv"
    temp_file="$OUTDIR/${chunk_name}.tmp"
    shared_log="$OUTDIR/mpileup.log"
    lockfile="$OUTDIR/mpileup.lock"

    # Skip if fully logged
    if grep -E -q "$chunk_name completed in [0-9]+s" "$shared_log" 2>/dev/null; then
        # echo "Skipping $chunk_name (already completed)"
        return
    fi

    # Remove any leftover temp or output
    rm -f "$output_file" "$temp_file"

    start_time=$(date +%s)

    # Run pipeline into temp file
    bcftools mpileup -f "$REF" -d 10000 -B -O u -R "$region_file" --annotate FORMAT/AD -b "$BAM" 2>/dev/null \
        | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' > "$temp_file"

    end_time=$(date +%s)
    duration=$((end_time - start_time))

    # Log completion first
    flock "$lockfile" bash -c "echo \"[$(date '+%Y-%m-%d %H:%M:%S')] $chunk_name completed in ${duration}s\" >> \"$shared_log\""

    # Rename temp to final output
    mv "$temp_file" "$output_file"
}

export -f run_mpileup
export REF BAM OUTDIR

# === HANDLE CTRL+C AND CLEANUP ===
cleanup_on_interrupt() {
    echo 'Interrupted. Cleaning up...'
    find "$OUTDIR" -name "*.tsv" | while read -r tsv_file; do
        chunk_name=$(basename "$tsv_file" .tsv)
        if ! grep -E -q "$chunk_name completed in [0-9]+s" "$OUTDIR/mpileup.log"; then
            echo "Removing incomplete file: $tsv_file"
            rm -f "$tsv_file"
        fi
    done
    exit 1
}

trap cleanup_on_interrupt SIGINT

# === CLEAN UP INCOMPLETE FILES BEFORE START ===
if [ -f "$OUTDIR/mpileup.log" ]; then
    find "$OUTDIR" -name "*.tsv" | while read -r tsv_file; do
        chunk_name=$(basename "$tsv_file" .tsv)
        if ! grep -E -q "$chunk_name completed in [0-9]+s" "$OUTDIR/mpileup.log"; then
            echo "Removing incomplete file: $tsv_file"
            rm -f "$tsv_file"
        fi
    done
fi

# === RUN IN PARALLEL ===
find "$CHUNK_DIR" -name "apos_exons_chr*_chunk_*.txt" | sort | parallel -j "$JOBS" run_mpileup
