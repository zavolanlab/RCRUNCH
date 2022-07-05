#!/bin/bash

# Tear down test environment
cleanup () {
    rc=$?
    rm -rf .cache/
    rm -rf .config/
    rm -rf .snakemake/
    rm -rf logs/
    rm -rf results/
    cd $user_dir
    echo "Exit status: $rc"
}
trap cleanup EXIT

# Set up test environment
set -eo pipefail  # ensures that script exits at first command that exits with non-zero status
set -u  # ensures that script exits when unset variables are used
# set -x  # facilitates debugging by printing out executed commands
user_dir=$PWD
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
cd $script_dir

if [[ -e "../input_files_test/Homo_sapiens.GRCh38.98.chr20.gtf" && -e "../input_files_test/Homo_sapiens.GRCh38.transcriptome.chr20.fa" && -e "../input_files_test/Homo_sapiens.GRCh38.dna_sm.chromosome.20.fa" ]]; then
    echo "input files already exist, 3/3"

else
    bash ../get_test_data.sh
    python ../get_test_subset.py
fi

snakemake  --unlock  --snakefile="../../workflow/Snakefile" --rerun-incomplete --cores 10 --local-cores 10 
# # Run tests
snakemake \
    --snakefile="../../workflow/Snakefile" \
    --profile="../../profiles/slurm-singularity" \
    --configfile="config.yaml" \
    --singularity-args="--bind $script_dir/../.." 

# Evaluate performance based on test results
echo "Evaluation of results"
python performance_evaluation.py

