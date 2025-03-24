# Script to run the sce2g snakemake pipeline

# Params:
# 1: Full path to cell_clusters.tsv file
# 2: Result output dir

snakemake \
--profile snakemake_slurm_profile \
--use-conda \
 --config \
 cell_clusters=$1 \
 results_dir=results/${2} \
 IGV_dir=results/$2
