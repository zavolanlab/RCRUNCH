which snakemake
which singularity
# /usr/bin/env
snakemake --unlock  --rerun-incomplete --cores 10 --local-cores 10 
# snakemake --dag -np --use-singularity | dot -Tpdf > dag.pdf
# snakemake --rulegraph --configfile config.yaml | dot -Tpng > rulegraph.png

snakemake \
-s workflow/Snakefile \
--profile profiles/slurm-singularity/ \
--configfile config.yaml