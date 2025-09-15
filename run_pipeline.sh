# Script to run the sce2g snakemake pipeline

# Params:
# 1: Full path to cell_clusters.tsv file
# 2: Result output dir

snakemake \
--profile snakemake_slurm_profile \
--use-conda \
--conda-frontend mamba \
--config \
 cell_clusters=$1 \
 results_dir=results/${2} \
 IGV_dir=results/$2 \
 gene_annotations="/maps/projects/cbmr_shared/people/wkq953/non-GDPR/segment/pipeline/resources/gencode.v32.annotation.gtf.gz" \
  make_IGV_tracks=True