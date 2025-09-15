RUN_TYPE_PROJECT="test.hWAT_all"

snakemake --rulegraph \
--profile snakemake_slurm_profile \
  --use-conda \
  --config \
  cell_clusters=/maps/projects/cbmr_shared/people/tqb695/GDPR/_targets/files/scE2G_cfg_clusters_w_peaks_file.all.${RUN_TYPE_PROJECT}.tsv \
  results_dir=results/rulegraph \
  IGV_dir=results/rulegraph \
  gene_annotations="resources/genome_annotations/gencode.v32.annotation.gtf.gz" \
  make_IGV_tracks=True | \
 dot -Tpng > rulegraph_${RUN_TYPE_PROJECT}.png
