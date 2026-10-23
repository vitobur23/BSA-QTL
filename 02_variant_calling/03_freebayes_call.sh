#!/bin/bash

# STEP 2 of 3: variant calling with FreeBayes on the two merged bulk BAMs.
#
# Each bulk BAM is treated as one pool (--pooled-continuous): reads are not forced into a
# small-integer genotype model. Only the raw read counts (RO, AO, DP) are used in the QTL
# analysis, and these do not depend on the declared ploidy, so no copy-number map or ploidy
# override is passed and ploidy is left at the FreeBayes default.
# The order of the -b options defines the sample order in the VCF: bsa_high, then bsa_low.
#
# FreeBayes on a whole genome is memory-intensive (we used a 2-core job with 256 GB RAM and
# about 60 hours of runtime).
#
# Requires: freebayes, samtools.
# Usage: bash 03_freebayes_call.sh

# --- EDIT THESE ----------------------------------------------------------
REF="/path/to/Arabidopsis_arenosa_genome.softmasked.fna"   # reference genome
BAMDIR="/path/to/bam"                                       # output of the mapping step
HIGH_BAM="${BAMDIR}/bsa_high.merged.bam"
LOW_BAM="${BAMDIR}/bsa_low.merged.bam"
OUT_DIR="/path/to/vcf"
OUT_VCF="bsa_arenosa_freebayes.raw.vcf"
# -------------------------------------------------------------------------

mkdir -p "$OUT_DIR"

# Index the reference and the BAM files if indexes are missing
[[ -f "${REF}.fai" ]]      || samtools faidx "$REF"
[[ -f "${HIGH_BAM}.bai" ]] || samtools index "$HIGH_BAM"
[[ -f "${LOW_BAM}.bai" ]]  || samtools index "$LOW_BAM"

echo "Running FreeBayes (pooled-continuous, default ploidy)..."
freebayes \
    --pooled-continuous \
    -n 4 \
    -f "$REF" \
    -b "$HIGH_BAM" \
    -b "$LOW_BAM" \
    -v "$OUT_DIR/$OUT_VCF"
FB_EXIT=$?

# Check the exit code and not only whether the file is non-empty: a run killed by the
# out-of-memory handler leaves a truncated VCF.
if [[ $FB_EXIT -eq 0 && -s "$OUT_DIR/$OUT_VCF" ]]; then
    echo "Done: $OUT_DIR/$OUT_VCF"
    echo "Variant count: $(grep -vc '^#' "$OUT_DIR/$OUT_VCF")"
else
    echo "FreeBayes did not finish successfully (exit code $FB_EXIT)."
    echo "The output is likely truncated. Check the log for an out-of-memory or kill message."
    exit 1
fi
