snakemake \
--profile snakemake_slurm_profile \
--use-conda \
--conda-frontend mamba \
--config \
 cell_clusters="/maps/projects/cbmr_shared/people/tqb695/GDPR/_targets/files/single_scE2G_cluster_file.real.muscle_15.tsv" \
 results_dir="results/muscle_test" \
 IGV_dir="results/muscle_test" \
 gene_annotations="/maps/projects/cbmr_shared/people/wkq953/non-GDPR/segment/pipeline/resources/gencode.v32.annotation.gtf.gz" \
 make_IGV_tracks=True