#!/bin/bash

# STEP 3 of 3: filter the raw FreeBayes VCF and index it.
#
# Filters:
#   - biallelic SNPs only
#   - mapping quality (MQM) >= 20
#   - strand-bias metrics SRP and SAP < 100
#   - per-bulk read depth (DP) between 4 and 1000
#   - per-bulk allele balance AO/(AO+RO) between 0.1 and 0.9
#
# Sample order in the VCF is bsa_high, then bsa_low (the -b order in the FreeBayes step),
# so sample index 0 = bsa_high and index 1 = bsa_low in the bcftools expressions below.
#
# Additional filters (minimum depth 4 in both bulks, depth difference between bulks,
# reference-allele frequency) are applied later in the QTLseqr script (04_qtl_analysis).
#
# Requires: bcftools, tabix.
# Usage: bash 04_filter_and_index.sh

# --- EDIT THESE ----------------------------------------------------------
WORKDIR="/path/to/vcf"                                  # folder with the raw VCF
RAW_VCF="${WORKDIR}/bsa_arenosa_freebayes.raw.vcf"
OUT_NAME="bsa_arenosa_freebayes_filtered"
# -------------------------------------------------------------------------

cd "$WORKDIR" || exit 1

# Keep biallelic SNPs only
bcftools view -v snps -m2 -M2 "$RAW_VCF" -Oz -o biallelic.vcf.gz
bcftools index -t biallelic.vcf.gz

# Apply quality, strand-bias, depth and allele-balance filters
bcftools view -i '
    INFO/MQM>=20
    && INFO/SRP<100 && INFO/SAP<100
    && FORMAT/DP[0]>=4 && FORMAT/DP[0]<=1000
    && FORMAT/DP[1]>=4 && FORMAT/DP[1]<=1000
    && (FORMAT/AO[0:0] / (FORMAT/AO[0:0] + FORMAT/RO[0])) >= 0.1
    && (FORMAT/AO[0:0] / (FORMAT/AO[0:0] + FORMAT/RO[0])) <= 0.9
    && (FORMAT/AO[1:0] / (FORMAT/AO[1:0] + FORMAT/RO[1])) >= 0.1
    && (FORMAT/AO[1:0] / (FORMAT/AO[1:0] + FORMAT/RO[1])) <= 0.9
' biallelic.vcf.gz -Oz -o "${OUT_NAME}.vcf.gz"
FILTER_EXIT=$?

if [[ $FILTER_EXIT -ne 0 ]]; then
    echo "bcftools filtering failed (exit code $FILTER_EXIT)."
    exit 1
fi

tabix -p vcf "${OUT_NAME}.vcf.gz"
TABIX_EXIT=$?

echo "Variant counts:"
echo -n "raw:        "; grep -vc "^#" "$RAW_VCF"
echo -n "biallelic:  "; bcftools view -H biallelic.vcf.gz | wc -l
FILTERED_COUNT=$(bcftools view -H "${OUT_NAME}.vcf.gz" 2>/dev/null | wc -l)
echo "filtered:   ${FILTERED_COUNT}"

if [[ $TABIX_EXIT -ne 0 || "$FILTERED_COUNT" -eq 0 ]]; then
    echo "Filtered VCF is empty or failed to index (tabix exit $TABIX_EXIT)."
    exit 1
fi

echo "Done: $WORKDIR/${OUT_NAME}.vcf.gz"
