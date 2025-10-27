#!/bin/bash
# Exit immediately if a command exits with a non-zero status.
set -e
# Treat unset variables as an error
set -u
# Exit on pipe failures
set -o pipefail

while getopts s:w:h:d:t: flag
do
    case "${flag}" in
        s) species=${OPTARG};;
        w) source=${OPTARG};;
        h) haploid=${OPTARG};;
        d) HOME=${OPTARG};;
        t) threads=${OPTARG};;
    esac
done

# --- Setup ---

# if no files are found, instead of expanding to the literal string "*.bam".
shopt -s nullglob

# Define directories
src_dir="${HOME}/${source}/${species}"
out_dir="${HOME}/aligned_bam/${species}"
final_bam="${out_dir}/${species}_merged_sorted.bam"

# Define reference
if [ "$haploid" == "True" ]; then
    ref="${HOME}/assemblies/${species}.pri.fasta"
else
    ref="${HOME}/assemblies/${species}.merged.fasta"
fi

echo "Source directory: $src_dir"
echo "Output directory: $out_dir"
echo "Reference: $ref"
echo "Final BAM: $final_bam"
echo "Threads: $threads"

# --- Checkpoint ---
# If the final bam file already exists, we are done.
if [ -f "$final_bam" ]; then
    echo "Final sorted BAM file already exists. Skipping."
    exit 0
fi

# --- Create Output Dirs ---
mkdir -p "${out_dir}"

# --- Find all read files ---
# We do this before streaming to know what we have
bam_files=(${src_dir}/*.bam)
fastq_gz_files=(${src_dir}/*.fastq.gz)
fastq_files=(${src_dir}/*.fastq)
fq_files=(${src_dir}/*.fq)

total_files=$(( ${#bam_files[@]} + ${#fastq_gz_files[@]} + ${#fastq_files[@]} + ${#fq_files[@]} ))

if [ $total_files -eq 0 ]; then
    echo "Error: No .bam, .fastq.gz, .fastq, or .fq files found in ${src_dir}" >&2
    exit 1
fi

echo "Found ${#bam_files[@]} BAM, ${#fastq_gz_files[@]} fastq.gz, ${#fastq_files[@]} fastq, and ${#fq_files[@]} fq files."

# --- Map and Sort ---
echo "--- Starting alignment stream ---"

(
    # Convert BAM files to FASTQ and stream to stdout
    if [ ${#bam_files[@]} -gt 0 ]; then
        echo "Streaming ${#bam_files[@]} BAM file(s)..."
        for f in "${bam_files[@]}"; do
            echo "Converting $f..."
            pbindex $f
            # -o /dev/stdout streams output instead of writing to a file
            bam2fastq -o /dev/stdout -j $threads $f
        done
    fi

    # Decompress and stream .fastq.gz files
    if [ ${#fastq_gz_files[@]} -gt 0 ]; then
        echo "Streaming ${#fastq_gz_files[@]} .fastq.gz file(s)..."
        zcat "${fastq_gz_files[@]}"
    fi

    # Stream plain .fastq files
    if [ ${#fastq_files[@]} -gt 0 ]; then
        echo "Streaming ${#fastq_files[@]} .fastq file(s)..."
        cat "${fastq_files[@]}"
    fi

    # Stream plain .fq files
    if [ ${#fq_files[@]} -gt 0 ]; then
        echo "Streaming ${#fq_files[@]} .fq file(s)..."
        cat "${fq_files[@]}"
    fi

) | minimap2 -t $threads -ax map-hifi $ref - | \
    samtools sort -@ $threads -o $final_bam -

echo "--- Mapping and sorting complete ---"
echo "Final file saved to: $final_bam"

# --- Index the sorted BAM file ---
echo "Indexing sorted BAM..."
samtools index -c -@ $threads $final_bam

echo "--- All done. ---"