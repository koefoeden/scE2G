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
 IGV_dir=results/$2 \
 gene_annotations="resources/genome_annotations/gencode.v32.annotation.gtf.gz" \
  make_IGV_tracks=True



snakemake \
--profile snakemake_slurm_profile \
--use-conda \
 --config \
 cell_clusters="/maps/projects/cbmr_shared/people/tqb695/GDPR/_targets/files/scE2G_cfg_clusters_w_peaks_file.all.real.muscle.tsv" \
 results_dir="results/2025-05-05-1017_muscle_colata_FLINC_PBMC" \
 IGV_dir="results/2025-05-05-1017_muscle_colata_FLINC_PBMC" \
 gene_annotations="resources/genome_annotations/gencode.v32.annotation.gtf.gz" \
 make_IGV_tracks=True
