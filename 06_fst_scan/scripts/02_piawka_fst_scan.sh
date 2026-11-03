#!/usr/bin/env bash
# Step 2: F_ST scan in 1-kb windows between the parental populations
# (non-serpentine VLA vs serpentine BOR) with piawka.
#
# Requires: bcftools, tabix, bgzip, and piawka (https://github.com/novikovalab/piawka,
#           v0.9.0, commit e524764). piawka must be on PATH.
#
# Inputs : filtered VCF from 01_filter_vcf.sh
#          input/groups_VLA_BOR.txt   (sample <TAB> population; 24 BOR + 23 VLA)
#          input/pop_pair.tsv         (population pair: VLA <TAB> BOR)
# Output : OUT/Fst_scan_<pop1>_<pop2>.bed.gz
#
# Usage: bash 02_piawka_fst_scan.sh arenosa_05_dp.snp.m.dm.vcf.gz output_dir

set -euo pipefail

VCF_IN="${1:?filtered VCF from 01_filter_vcf.sh}"
OUT="${2:?output directory}"
GROUPS_FILE="input/groups_VLA_BOR.txt"
PAIRS_FILE="input/pop_pair.tsv"

WINDOW=1000          # window size (bp)
MISS=0.5             # sites with more than 50 % missing genotypes in a population are removed
LINES=10             # windows with fewer than 10 SNPs (lines) in a population are removed
JOBS=8

mkdir -p "$OUT"

## non-overlapping 1-kb windows; scaffold sizes are read from the VCF header
piawka win -v "$VCF_IN" -w "$WINDOW" > "$OUT/windows_1kb.bed"

while IFS=$'\t' read -r pop1 pop2 || [[ -n "$pop1" ]]; do
    echo "Processing $pop1 + $pop2"

    ## samples of the pair; keep only variable sites, which makes the SNP count per window correct
    grep -P "\t$pop1$|\t$pop2$" "$GROUPS_FILE" | cut -f1 > "$OUT/samples_${pop1}_${pop2}.txt"
    PAIR_VCF="$OUT/pair_${pop1}_${pop2}.vcf.gz"
    bcftools view -S "$OUT/samples_${pop1}_${pop2}.txt" -Ou "$VCF_IN" |
        bcftools +fill-AN-AC -Ou |
        bcftools view -a -Ou |
        bcftools view -m2 -o "$PAIR_VCF"
    bcftools index --tbi "$PAIR_VCF"

    ## per-site statistics -> filter sites -> sum within windows -> filter windows
    piawka calc --persite --jobs $JOBS --vcf "$PAIR_VCF" \
            --groups <(grep -P "\t$pop1$|\t$pop2$" "$GROUPS_FILE") --bed "$OUT/windows_1kb.bed" \
            --stats fst,fstwc,afd,miss,lines,segr --mult --dependencies --quiet |
        piawka filt -e "miss <= $MISS" |
        piawka sum --bed "$OUT/windows_1kb.bed" |
        piawka filt -e "lines >= $LINES" |
        bgzip > "$OUT/Fst_scan_${pop1}_${pop2}.bed.gz"
done < "$PAIRS_FILE"

echo "Done."
