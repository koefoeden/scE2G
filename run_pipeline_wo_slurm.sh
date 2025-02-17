source ~/.bashrc
init_mamba
conda activate sc-E2G

cd ~/tqb695/sc-E2G

 snakemake \
 --config cell_clusters="/home/tqb695/tqb695/sc-E2G/config/config_cell_clusters_pipeline.tsv" \
 results_dir="/home/tqb695/tqb695/sc-E2G/results/pipeline/" \
 IGV_dir="/home/tqb695/tqb695/sc-E2G/results/pipeline/" \
--default-resources runtime=3600 \
 --jobs 2 \
 --cores 10 \
 --use-conda \
 --rerun-incomplete