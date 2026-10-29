# BSA-QTL
We combined bulk segregant analysis (BSA) of an F2 mapping population (serpentine x non-serpentine A. arenosa) with low-coverage individual sequencing, and tested the resulting QTL regions against candidate genes from a genome-wide FST divergence scan of the parental populations.

Article: [DOI TO BE ADDED]
Raw sequencing data: European Nucleotide Archive, accession PRJEB111765
Reference genome and annotation: https://github.com/PaulGrini/Arabidopsis_assemblies/tree/main/01_arenosa_assembly

## Analysis overview


| 1 | Read mapping to the *A. arenosa* reference (Snakemake) 
| 2 | Pooled variant calling (FreeBayes) on the two bulk alignments, filtering 
| 3 | QTL analysis (QTLseqr, G' and delta-SNP), gene assignment, overlap with FST candidates 
| 4 | GO enrichment of the overlapping genes (topGO) 
| 5 | Heritability of rosette size (full-sib design) 
| 6 | FST divergence scan of the parental populations 
 




Steps are numbered in the order in which they were run.

## Requirements

- Snakemake and the workflow from <https://github.com/jgerchen/polyploid_variant_calling>
  (commit f794f85) for steps 1 and 8
- FastQC v0.12.1, Trimmomatic v0.39, BWA-MEM v0.7.18, Picard v3.1.1
- FreeBayes, bcftools, samtools, tabix (versions: [TO BE FILLED IN])
- R [VERSION] with QTLseqr v0.7.5.2, topGO v2.64.0, SuperExactTest, data.table, dplyr, tidyr,
  readr, ggplot2, cowplot, zoo, rlang


## Step details

### 1. Read mapping (`01_mapping/`)
Follow workflow from <https://github.com/jgerchen/polyploid_variant_calling>
`bsa.yaml` and `bsa_samples.tsv` define the Snakemake run. Each of the two rows
(`high`, `low`) lists the sequencing codes of all 62 individuals of that bulk. Reads are trimmed,
mapped to the soft-masked *A. arenosa* assembly, duplicate-marked, and the individual alignments are
merged into one BAM file per bulk. Run with `run_snake_bsa_1_mapping.sh`.

### 2. Variant calling (`02_variant_calling/`)
Run in order:

1. `01_verify_bams.sh` checks the merged bulk BAMs.
2. `03_freebayes_call.sh` calls variants jointly on the two bulk BAMs with
   `--pooled-continuous`. Only the read-count fields (RO, AO, DP) are used downstream;
   genotype calls and ploidy-derived fields are not used.
3. `04_filter_and_index.sh` keeps biallelic SNPs with mapping quality >= 20, SRP and SAP < 100, and
   per-bulk depth between 4 and 1000.


Output: `bsa_arenosa_freebayes_filtered.vcf.gz` (2,161,257 SNPs; too large for this repository)

### 3. QTL analysis (`03_qtl_analysis/`)
`QTLseqr_BSA_final_ws1000000_DP4_top1pct.R` reads the filtered VCF and:

- keeps SNPs with depth >= 4 in both bulks, depth difference <= 100 and reference-allele
  frequency between 0.05 and 0.95
- computes SNP indices and delta-SNP from read counts
- runs `runQTLseqAnalysis()` (F2, 1 Mb window, 248 chromosomes per bulk, 10,000 replications) and
  `runGprimeAnalysis()` (1 Mb window), FDR 0.05
- assigns SNPs to genes (gene body +/- 2 kb, clipped at the neighbouring gene)
- tests the overlap with FST candidate genes (Fisher's exact test, SuperExactTest)

The script expects the filtered VCF and the FST candidate table (`input/`) at the paths set at its
top. Output files are tagged `ws1000000_DP4_arenosa_top1pct`.

### 4. GO enrichment (`04_go_enrichment/`)
`GO_QTLseqr_arenosa_rerun_top1pct.R` runs topGO (Fisher's exact test, `elim` algorithm) on the genes
in the BSA x FST overlap, against all annotated genes on the eight assembled scaffolds. GO terms
are transferred from *A. thaliana* orthologs (`GO_mappings_Arabidopsis_arenosa.tsv`). The
enrichment is run on the overlap set only.

### 5. Heritability (`05_heritability/`)
`heritability_fullsib_ANOVA.R` groups F2 individuals into full-sib families by their (unordered)
F1 parental pair, partitions the variance with a one-way ANOVA, and reports the full-sib intraclass
correlation t and the bound h2 <= 2t with family-level bootstrap confidence intervals. It also
tests the reciprocal-cross effect. Inputs: `crossing_design.csv` / genealogy table and the
phenotype table (`F2-both-treatments.xlsx`).

### 6. FST scan (`06_fst_scan/`)
[TO BE ADDED: FST in 1 kbp windows with piawka, top 1% windows, window-to-gene assignment with
2 kbp flanks.] Genotype calling for the parental populations follows the workflow cited above.

### 7. Tissue-specific expression (07_gene_expression/)

Gene_expression_BSA-NEW.R classifies the candidate genes by expression level in root tip, root and leaf tissues (rosette, cauline and adult vascular leaves) using the Arabidopsis tissue atlas (Expression Atlas E-MTAB-7978). Genes of the whole atlas are classified in each tissue as lowly (<= 25th percentile of mean TPM), medium or highly (>= 75th percentile) expressed, and the proportions of the candidate genes in these classes are compared with the genome-wide proportions expected for a random set of genes.

Inputs (put into input/):

E-MTAB-7978-query-results.tpms.tsv: TPM table, download from https://www.ebi.ac.uk/gxa/experiments/E-MTAB-7978/Downloads
Candidate_genes_AT_codes.txt: tab-separated table with a column A.thaliana containing the A. thaliana orthologs of the candidate genes (from Supplementary Data S2)

## Input data not in this repository

- Raw reads: ENA PRJEB111765 , [[[[FST data still to add]]]]
- Reference genome and annotation: see link above
- A. thaliana ortholog table: <https://github.com/mduchoslav/Brassicaceae_orthology>
- Phenotype data: see Supplementary Data of the article

## Notes on paths

Scripts were written for a local layout and use relative paths from the repository root. If a script
fails on a missing input, check the path variables at its top.


## Contact

Veronika Lipánová, lipan.veronika@gmail.com
