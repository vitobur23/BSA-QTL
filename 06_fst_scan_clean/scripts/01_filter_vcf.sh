#!/usr/bin/env bash
# Step 1: Filter the joint A. arenosa VCF used for the F_ST divergence scan.
#
# Input: joint VCF of all samples called with the polyploid_variant_calling
#        Snakemake workflow (https://github.com/jgerchen/polyploid_variant_calling)
#        against the A. arenosa reference genome, see config/variant_calling.yaml.
# Requires: bcftools, tabix, bedtools, python3, and make_depth_mask.py from
#        https://github.com/jgerchen/polyploid_popgen/tree/main/depth_mask
#
# Filtering steps (named after the intermediate files):
#   01  tag sites that were originally multiallelic
#   02  genotypes with FMT/DP <= MIN_DP set to missing; unobserved alleles and
#       invariant sites removed
#   03  keep SNPs only (no indels, MNPs or sites overlapping deletions)
#       + keep only the A. arenosa samples listed in SAMPLE_LIST
#   04  remove sites with more than MAX_MISSING of genotypes missing
#   05  remove sites in the depth mask (excess depth, e.g. collapsed repeats)
#
# Usage: bash 01_filter_vcf.sh joint.vcf.gz sample_list.txt make_depth_mask.py output_dir

set -euo pipefail

VCF_IN="${1:?joint VCF (bgzipped, indexed)}"
SAMPLE_LIST="${2:?file with one sample ID per line (A. arenosa samples retained)}"
DEPTH_MASK_SCRIPT="${3:?path to make_depth_mask.py}"
OUT="${4:?output directory}"

MIN_DP=3                 # genotypes with DP <= 3 are set to missing (DP >= 4 retained)
MAX_MISSING=0.5          # maximal fraction of missing genotypes per site
MASK_PROP_EXCESS=0.5     # fraction of samples with excess depth needed to mask a site
MASK_MERGE_DISTANCE=500  # bp; nearby masked sites are merged into regions
THREADS=2

mkdir -p "$OUT"
BASE="arenosa"

## 01: tag originally multiallelic sites
bcftools filter --threads $THREADS -m+ -s 'MultiAllelic' -e 'N_ALT>1' \
    -o "$OUT/${BASE}_01_tagMultiAllelic.vcf.gz" "$VCF_IN"

## 02: low-depth genotypes -> missing, drop unobserved alleles and invariant sites
bcftools view -h "$OUT/${BASE}_01_tagMultiAllelic.vcf.gz" > "$OUT/${BASE}_01.header"
sed -i 's/ID=PL,Number=G/ID=PL,Number=./g' "$OUT/${BASE}_01.header"   # header fix required by the next command
bcftools reheader -h "$OUT/${BASE}_01.header" "$OUT/${BASE}_01_tagMultiAllelic.vcf.gz" |
    bcftools filter -i "FMT/DP>$MIN_DP" --set-GTs . --threads $THREADS -Ou |
    bcftools +fill-AN-AC -Ou |
    bcftools view -a -Ou |
    bcftools view -m2 --threads $THREADS -o "$OUT/${BASE}_02_dp.vcf.gz"
tabix "$OUT/${BASE}_02_dp.vcf.gz"

## 03: SNPs only; keep the selected samples; drop alleles absent from them and invariants
bcftools view -m2 -e 'ALT="*"' --exclude-types indels,mnps --threads $THREADS -Ou "$OUT/${BASE}_02_dp.vcf.gz" |
    bcftools view --threads $THREADS -S "$SAMPLE_LIST" -Ou |
    bcftools +fill-AN-AC -Ou |
    bcftools view -a -Ou |
    bcftools view -m2 --threads $THREADS -o "$OUT/${BASE}_03_dp.snp.vcf.gz"
tabix "$OUT/${BASE}_03_dp.snp.vcf.gz"

## 04: site missingness
bcftools view -i "F_MISSING<=$MAX_MISSING" --threads $THREADS \
    -o "$OUT/${BASE}_04_dp.snp.m.vcf.gz" "$OUT/${BASE}_03_dp.snp.vcf.gz"
tabix "$OUT/${BASE}_04_dp.snp.m.vcf.gz"

## depth mask: sites where at least MASK_PROP_EXCESS of samples have excess depth
N_LOCI=$(zcat "$OUT/${BASE}_04_dp.snp.m.vcf.gz" | grep -c '^[^#]')
N_SAMPLES=$(bcftools query -l "$OUT/${BASE}_04_dp.snp.m.vcf.gz" | wc -l)
N_EXCESS=$(echo "$N_SAMPLES * $MASK_PROP_EXCESS" | bc)
mkdir -p "$OUT/depthmask"
python3 "$DEPTH_MASK_SCRIPT" \
    -v "$OUT/${BASE}_04_dp.snp.m.vcf.gz" \
    -o "$OUT/depthmask/${BASE}.depthmask.tsv" \
    -c "$OUT/depthmask/${BASE}.depthmask_counts.tsv" \
    -p "$OUT/depthmask/${BASE}.depthmask_hist.tsv" \
    -n "$N_EXCESS" -l "$N_LOCI" -m "$MIN_DP" > "$OUT/depthmask/${BASE}.depthmask_script_output.txt"
awk '{printf "%s\t%d\t%d\n", $1, ($2 - 1), $2}' "$OUT/depthmask/${BASE}.depthmask.tsv" > "$OUT/depthmask/${BASE}.depthmask.bed"
bedtools merge -d "$MASK_MERGE_DISTANCE" -i "$OUT/depthmask/${BASE}.depthmask.bed" > "$OUT/depthmask/${BASE}.depthmask_merged.bed"

## 05: apply the depth mask
bcftools view --threads $THREADS -T "^$OUT/depthmask/${BASE}.depthmask_merged.bed" \
    -o "$OUT/${BASE}_05_dp.snp.m.dm.vcf.gz" "$OUT/${BASE}_04_dp.snp.m.vcf.gz"
tabix "$OUT/${BASE}_05_dp.snp.m.dm.vcf.gz"

echo "Final VCF for the F_ST scan: $OUT/${BASE}_05_dp.snp.m.dm.vcf.gz"
