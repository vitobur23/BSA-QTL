# 06_fst_scan: F_ST divergence scan of the parental populations

Genome-wide scan of differentiation between the non-serpentine (VLA, N1) and serpentine
(BOR, S1) parental populations of *A. arenosa* in 1-kb windows, mapped to the native
*A. arenosa* reference genome. The top 1 % of windows are assigned to genes and
define the F_ST candidate genes that are overlapped with the BSA-QTL regions
(`../03_qtl_analysis`).

Samples: 24 individuals from BOR (S) and 23 from VLA (N), see `input/groups_VLA_BOR.txt`
(8 + 8 individuals from Konecna et al. 2021, 1 + 1 from Celestini et al., and the newly
sequenced individuals; accession codes in the Data availability statement of the paper).

## Workflow

| Step | File | What it does |
|---|---|---|
| 0 | `config/variant_calling.yaml` | Settings of the Snakemake variant-calling workflow ([polyploid_variant_calling](https://github.com/jgerchen/polyploid_variant_calling); reads of all *A. arenosa* samples mapped to the *A. arenosa* genome, GATK joint calling). Run as described in that repository. |
| 1 | `scripts/01_filter_vcf.sh` | Filters the joint VCF: genotypes with DP < 4 set to missing, SNPs only, sites with > 50 % missing genotypes removed, depth mask applied. |
| 2 | `scripts/02_piawka_fst_scan.sh` | F_ST (Weir & Cockerham estimator, `fstwc`, and the plain `fst`) in 1-kb windows with [piawka](https://github.com/novikovalab/piawka) v0.9.0 (commit e524764): sites with > 50 % missing genotypes in a population are dropped, windows need >= 10 SNPs in both populations. |
| 3 | `scripts/03_fst_candidates.R` | Ranks windows by F_ST (`fstwc`), assigns the windows to genes (gene body + 2 kb flanks), writes `fst_candidates.tsv` and the list of unique candidate genes. |

## Inputs

- Reference genome: `Arabidopsis_arenosa_genome.softmasked.fna` and annotation `Arabidopsis_arenosa_genome.annotation.gtf`
  (Bramsiepe et al. 2023, <https://github.com/PaulGrini/Arabidopsis_assemblies>).
- Joint VCF of all *A. arenosa* samples (output of step 0) and a list of the samples retained in step 1.
- `input/groups_VLA_BOR.txt`, `input/pop_pair.tsv` (provided).

## Output

`results/fst_candidates.tsv` is provided in this repository (output of step 3). It holds one row per
1-kb window x gene region (gene body, upstream region, downstream region) for windows above the 95th
percentile of `fstwc`. Windows with no gene have `NA` in `gene_ID`. Only rows with
`fstwc.percentile >= 99` (top 1 % of windows) and scaffolds 1-8 are used as F_ST candidates: this yields the
1,162 candidate genes (`../03_qtl_analysis` applies this filter).

## Notes

- Percentiles are computed over all windows that pass the SNP-number and missingness filters.
  Because a window can overlap several gene regions, it appears in several rows of `fst_candidates.tsv`;
  the number of rows is therefore larger than the number of unique windows.
- Flanking regions are measured from the start and stop codon (the annotation has no UTRs) and are
  shortened where they would reach the neighbouring gene.
