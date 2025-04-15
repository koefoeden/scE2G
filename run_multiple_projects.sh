#!/bin/bash
# Usage: ./script.sh project1 project2 ...

# Check that at least one project name was provided.
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 project1 project2 ..."
    exit 1
fi

# Read project names from command-line arguments.
RUN_TYPE_PROJECTS=("$@")

# Concatenate project names with underscores to form a unique identifier.
RESULT_NAME=results/$(date +%Y-%m-%d-%H%M)_$(IFS=_; echo "${RUN_TYPE_PROJECTS[*]}")

# Define the output concatenated cell_clusters file.
OUTPUT_CELL_CLUSTERS=${RESULT_NAME}.tsv


# Concatenate TSV files:
# Assumes each TSV file has the same header; we keep the header from the first file only.
FIRST=1
for proj in "${RUN_TYPE_PROJECTS[@]}"; do
  INPUT_FILE="/maps/projects/cbmr_shared/people/tqb695/GDPR/_targets/files/scE2G_cfg_clusters_w_peaks_file.all.real.${proj}.tsv"
  if [ $FIRST -eq 1 ]; then
    # Write the header from the first file.
    head -n 1 "$INPUT_FILE" > "$OUTPUT_CELL_CLUSTERS"
    # Process data lines: prepend project name to the first column.
    tail -n +2 "$INPUT_FILE" | awk -v proj="${proj}" 'BEGIN { FS=OFS="\t" } { $1 = proj "_" $1; print }' >> "$OUTPUT_CELL_CLUSTERS"
    FIRST=0
  else
    # For subsequent files, skip the header and process data lines.
    tail -n +2 "$INPUT_FILE" | awk -v proj="${proj}" 'BEGIN { FS=OFS="\t" } { $1 = proj "_" $1; print }' >> "$OUTPUT_CELL_CLUSTERS"
  fi
done

# exit 0
# Run snakemake once using the concatenated cell_clusters file and combined results directory.

snakemake \
  --profile snakemake_slurm_profile \
  --use-conda \
  --config \
  cell_clusters="$OUTPUT_CELL_CLUSTERS" \
  results_dir=${RESULT_NAME} \
  IGV_dir=${RESULT_NAME} \
  gene_annotations="resources/genome_annotations/gencode.v32.annotation.gtf.gz" \
  make_IGV_tracks=True
