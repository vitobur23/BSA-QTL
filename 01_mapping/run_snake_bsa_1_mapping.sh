#!/bin/bash

### Activate the conda environment that contains Snakemake (here named "Snakemake")
source /path/to/mambaforge/bin/activate Snakemake

### Change this to the folder where your workflow is located
cd /path/to/polyploid_variant_calling/workflow

# This is the snakemake command for read mapping. It will map reads for all individuals in the samples list and merge the alignments of each bulk.
# The option --rerun-triggers mtime ensures that Snakemake will only generate output files that do not exist yet, since it likes to rerun rules where the code has been modified.
# Change --cores to the number of CPU cores you want to use.
# On a cluster, replace "--cores 8" with your Snakemake profile, e.g. "--profile /path/to/snakemake_profile -j 100 --max-status-checks-per-second 1".

snakemake /path/to/output/report/merge_bam_stats/bsa.bam.table.tsv --configfile /path/to/bsa.yaml --cores 8 --rerun-triggers mtime