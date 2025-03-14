# Script to run the sce2g snakemake pipeline

PROJECT=$1

snakemake --dag all \
--profile snakemake_slurm_profile \
--use-conda \
 --config \
 cell_clusters=/maps/projects/cbmr_shared/people/tqb695/GDPR/_targets/objects/scE2G_cfg_clusters_w_peaks_file.subAll.real.${PROJECT}.tsv \
 results_dir=results/${PROJECT} \
 IGV_dir=results/${PROJECT} | dot -Tpdf > ${PROJECT}_dag.pdf



snakemake --rulegraph all \
--profile snakemake_slurm_profile \
--use-conda \
 --config \
 cell_clusters=/maps/projects/cbmr_shared/people/tqb695/GDPR/_targets/objects/scE2G_cfg_clusters_w_peaks_file.subAll.real.${PROJECT}.tsv \
 results_dir=results/${PROJECT} \
 IGV_dir=results/${PROJECT} | dot -Tpdf > ${PROJECT}_rulegraph.pdf
