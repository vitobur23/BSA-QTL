#!/bin/bash

# STEP 1 of 3: verify the merged bulk BAMs before variant calling.
#
# We check that each merged BAM contains all individuals of its bulk (one library, LB,
# per individual) and that all read groups share the same sample name (SM). With a
# single SM per file, FreeBayes treats each BAM as ONE pooled sample ("bsa_high" /
# "bsa_low"), so no read-group rewriting is needed.
#
# Requires: samtools (and GNU grep for the -P option).
# Usage: bash 01_verify_bams.sh

# --- EDIT THESE ----------------------------------------------------------
BAMDIR="/path/to/bam"                       # output of the mapping step (01_mapping)
HIGH_BAM="${BAMDIR}/bsa_high.merged.bam"
LOW_BAM="${BAMDIR}/bsa_low.merged.bam"
EXPECTED_INDIVIDUALS=62                     # individuals per bulk
# -------------------------------------------------------------------------

check_bam () {
    local label="$1" bam="$2"
    echo "=== $label: $bam ==="

    local sm_values
    sm_values=$(samtools view -H "$bam" | grep "^@RG" | grep -oP 'SM:\K[^\t]+' | sort -u)
    local sm_count
    sm_count=$(echo "$sm_values" | grep -c .)
    echo "Distinct SM (sample) values: $sm_count -> $sm_values"
    if [[ "$sm_count" -ne 1 ]]; then
        echo "WARNING: expected exactly one SM value per bulk BAM -- FreeBayes would"
        echo "         otherwise see this file as more than one sample."
    fi

    local lb_count
    lb_count=$(samtools view -H "$bam" | grep "^@RG" | grep -oP 'LB:\K[^\t]+' | sort -u | wc -l)
    echo "Distinct LB (individual) values: $lb_count (expected $EXPECTED_INDIVIDUALS)"
    if [[ "$lb_count" -ne "$EXPECTED_INDIVIDUALS" ]]; then
        echo "WARNING: individual count does not match bsa_samples.tsv -- check before"
        echo "         running FreeBayes."
    fi
    echo
}

check_bam "HIGH bulk" "$HIGH_BAM"
check_bam "LOW bulk"  "$LOW_BAM"
