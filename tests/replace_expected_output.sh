#!/bin/bash
# Used when we want to change the expected output to match the latest test output run
set -e

if [ $# -ne 1 ]; then
    echo "Usage: $0 <results_dir>"
    exit 1
fi
results_dir=$1

# remove old expected output
rm -rf expected_output/powerlaw_v2_models/*

# replace with new files except full predictions (too big)
mkdir -p powerlaw_v2_models/K562_chr22_cluster1/multiome_powerlaw_v2/
cp -R "${results_dir}/K562_chr22_cluster1/multiome_powerlaw_v2/"* powerlaw_v2_models/K562_chr22_cluster1/multiome_powerlaw_v2/
rm  "${results_dir}/K562_chr22_cluster1/multiome_powerlaw_v2/"* powerlaw_v2_models/K562_chr22_cluster1/multiome_powerlaw_v2/encode_re2g_predictions.tsv.gz