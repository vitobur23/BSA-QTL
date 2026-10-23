# QTL analysis of the BSA data with QTLseqr (G' and delta-SNP approaches), assignment of
# significant SNPs to genes, and overlap of the QTL genes with FST candidate genes.
#
# Run from the folder that contains this script (03_qtl_analysis) after putting the input
# files into ./input (see README):
#   input/bsa_arenosa_freebayes_filtered.vcf.gz       filtered VCF from 02_variant_calling
#   input/Arabidopsis_arenosa_genome.annotation.gtf.gz gene annotation of the A. arenosa genome
#   input/Arabidopsis_arenosa_At_orthologues_and_GO.tsv A. arenosa - A. thaliana ortholog table
#   input/fst_candidates.tsv                           FST divergence scan windows (Supplementary Data S3)
# All output files are written to the working directory.
#setwd("/path/to/this/folder")

# Install required packages (run once)
# install.packages(c("data.table","dplyr","readr","tidyr","cowplot","zoo","rlang",
#                     "ape","ggvenn"))
#devtools::install_github("bmansfeld/QTLseqr")
# BiocManager::install("SuperExactTest")

library(QTLseqr)
library(data.table)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(zoo)
library(rlang)

##### Initial Data Preparation #####

# Define bulk names (labels used only for plot legends)
highBulk <- "highbulk"
lowBulk  <- "lowbulk"

# Every output file from this script is tagged with this string, so window size
# and depth threshold are always visible in the file names.
tag <- "ws1000000_DP4_arenosa_top1pct"

# The VCF is the FreeBayes joint call of the two merged bulk BAMs (--pooled-continuous),
# restricted to biallelic SNPs and filtered with bcftools in 02_variant_calling (see the
# ##bcftools_view* header lines in the VCF for the exact filters).
# Samples in the VCF are named "bsa_high" and "bsa_low" (in that column order).
vcf_file <- "input/bsa_arenosa_freebayes_filtered.vcf.gz"

# Read directly from the VCF instead of a GATK VariantsToTable export (that
# intermediate step is no longer needed now that we work straight off the FreeBayes
# VCF). Requires `gzip` on PATH to decompress on the fly (present via Rtools/Git on
# Windows, native on macOS/Linux); if that's not available, install "R.utils" and
# call fread(vcf_file, ...) directly instead of via cmd=.
vcf <- fread(cmd = paste("gzip -dc", shQuote(vcf_file)), skip = "#CHROM", header = TRUE,
             select = c("#CHROM", "POS", "REF", "ALT", "FORMAT", "bsa_high", "bsa_low"))
setnames(vcf, "#CHROM", "CHROM")

# FORMAT is GT:DP:RO:QR:AO:QA:GL for every record in this (--pooled-continuous) VCF
# (verified on the full file) - i.e. DP is field 2, RO field 3, AO field 5. Guard
# against that assumption silently breaking if the VCF is ever regenerated with a
# different field order.
stopifnot(all(vcf$FORMAT == "GT:DP:RO:QR:AO:QA:GL"))

to_num <- function(x) { x[x == "."] <- NA; as.numeric(x) }

tmp_high <- tstrsplit(vcf$bsa_high, ":", fixed = TRUE)
tmp_low  <- tstrsplit(vcf$bsa_low,  ":", fixed = TRUE)

vcf[, `:=`(
  DP_highbulk = to_num(tmp_high[[2]]),
  RO_highbulk = to_num(tmp_high[[3]]),
  AO_highbulk = to_num(tmp_high[[5]]),
  DP_lowbulk  = to_num(tmp_low[[2]]),
  RO_lowbulk  = to_num(tmp_low[[3]]),
  AO_lowbulk  = to_num(tmp_low[[5]])
)]

# FreeBayes doesn't emit GQ (Genotype Quality) even though the VCF header declares
# the FORMAT field - it's simply never populated per-record - so GQ is set to a
# constant "good" value below.
file <- as.data.frame(vcf[, .(CHROM, POS, REF, ALT,
                               DP_lowbulk, RO_lowbulk, AO_lowbulk,
                               DP_highbulk, RO_highbulk, AO_highbulk)])
rm(vcf, tmp_high, tmp_low)

##### Coverage Histograms #####

ggplot(file, aes(x = DP_lowbulk)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  theme_minimal() +
  labs(title = "DP Histogram - lowbulk", x = "Depth (DP_lowbulk)", y = "Variant count")

ggplot(file, aes(x = DP_highbulk)) +
  geom_histogram(binwidth = 5, fill = "tomato", color = "black") +
  theme_minimal() +
  labs(title = "DP Histogram - highbulk", x = "Depth (DP_highbulk)", y = "Variant count")

##### Filtering #####
# Minimum depth of 4 in both bulks, matching the floor of the upstream bcftools filter
# (FORMAT/DP>=4, see 02_variant_calling/04_filter_and_index.sh).

file_dp4 <- file %>%
  filter(DP_lowbulk >= 4, DP_highbulk >= 4) %>%
  filter(abs(DP_highbulk - DP_lowbulk) <= 100) # added according to Mansfeld and Grumet 2018 and Mansfeld 2019 tutorial

file_dp4 <- file_dp4 %>%
  mutate(REF_FRQ = (RO_highbulk + RO_lowbulk)/(DP_lowbulk + DP_highbulk)) %>%
  filter(REF_FRQ <= 0.95, REF_FRQ >= 0.05) # added to calculate ref allele frequency and apply the filter according to Mansfeld and Grumet 2018

# Adding GQ values, Genotype quality - freebayes is not creating this data so I need to set it here, that all of the quality is good
file_dp4$GQ_lowbulk <- 99
file_dp4$GQ_highbulk <- 99

##### Manual SNP Index Calculation #####
# The QTLseqr import functions are written for GATK output, so the SNP indices and delta-SNP are calculated manually from the FreeBayes read counts (RO, AO, DP).

df_4 <- file_dp4 %>%
  filter(!is.na(DP_lowbulk), !is.na(AO_lowbulk),
         !is.na(DP_highbulk), !is.na(AO_highbulk)) %>%
  mutate(
    SNPindex_lowbulk = AO_lowbulk / DP_lowbulk,
    SNPindex_highbulk = AO_highbulk / DP_highbulk,
    deltaSNP = SNPindex_highbulk - SNPindex_lowbulk,
    DP_diff = abs(DP_highbulk - DP_lowbulk)
  )

Chroms <- paste0("scaffold_", 1:8)
df_4 <- df_4 %>% filter(CHROM %in% Chroms)

# Rename for compatibility

df_4 <- df_4 %>%
  mutate(
    AD_REF.LOW = RO_lowbulk,
    AD_ALT.LOW = AO_lowbulk,
    AD_REF.HIGH = RO_highbulk,
    AD_ALT.HIGH = AO_highbulk,
    DP.LOW = DP_lowbulk,
    DP.HIGH = DP_highbulk
  )

# check distributions after filtrations
summary(df_4$DP_lowbulk)
hist(df_4$DP_lowbulk)
summary(df_4$DP_highbulk)
hist(df_4$DP_highbulk)

##### QTLseq and G' Analysis #####

df_filt <- runQTLseqAnalysis(
  SNPset = df_4,
  windowSize = 1000000,
  popStruc = "F2", # according to crossings
  bulkSize = c(248, 248), # 62 tetraploid individuals = 248 chromosome copies - independent of the FreeBayes call itself (--pooled-continuous, no ploidy/cnv-map), this only shapes the QTLseqr null-distribution simulation
  replications = 10000,
  intervals = c(95, 99)
)

df_filt <- runGprimeAnalysis(
  SNPset = df_filt,
  windowSize = 1000000,
  outlierFilter = "deltaSNP"
)

##### Plotting #####

pdf(paste0("QTLseqr_plots_", tag, ".pdf"), width = 10, height = 7)

plotGprimeDist(SNPset = df_filt, outlierFilter = "Hampel")
plotGprimeDist(SNPset = df_filt, outlierFilter = "deltaSNP", filterThreshold = 0.05)

plotQTLStats(SNPset = df_filt, var = "Gprime", plotThreshold = TRUE, q = 0.05)
plotQTLStats(SNPset = df_filt, var = "deltaSNP", plotIntervals = TRUE)

QTLplots <- plotQTLStats(
  SNPset = df_filt,
  var = "negLog10Pval",
  plotThreshold = TRUE,
  q = 0.05
)
QTLplots

dev.off()

##### Extract Results #####

QTL_regions <- getSigRegions(SNPset = df_filt, alpha = 0.05)
all_sig_snps <- bind_rows(QTL_regions)

# Diagnostics: FDR threshold, smallest p-value and per-chromosome maxima, used to check
# which chromosomes carry significant windows.
fdrT <- QTLseqr::getFDRThreshold(df_filt$pvalue, alpha = 0.05)
fdrT
min(df_filt$pvalue, na.rm = TRUE)
-log10(min(df_filt$pvalue, na.rm = TRUE))

df_filt %>%
  group_by(CHROM) %>%
  summarise(maxGprime = max(Gprime, na.rm = TRUE),
            minPval   = min(pvalue, na.rm = TRUE),
            maxAbsDeltaSNP = max(abs(deltaSNP), na.rm = TRUE)) %>%
  arrange(minPval)

summary(df_filt$Gprime)

write_csv(all_sig_snps, paste0("BSA_QTLseqr_significant_SNPs_", tag, "_FDR0.05.csv"))

qtl_summary <- getQTLTable(SNPset = df_filt, method = "Gprime", alpha = 0.05)

getQTLTable(
  SNPset = df_filt,
  alpha = 0.05,
  export = TRUE,
  fileName = paste0("BSA_QTLanalysis_significant_QTLs_", tag, ".csv")
)

# get close up graph for the significant chromosomes - derived from the actual
# significant SNPs found above instead of a hardcoded list of scaffolds.
sig_chroms <- sort(unique(all_sig_snps$CHROM))
QTLplots <- plotQTLStats(
  SNPset = df_filt,
  var = "negLog10Pval",
  plotThreshold = TRUE,
  q = 0.05,
  subset = sig_chroms
)
QTLplots

############# Additional graphs ################################################
# plotting delta SNP (absolutni hodnota) and SNP indexes in significant QTL regions

# Helper to find the correct column names in the QTL table
std_col <- function(df, ...) {
  cands <- rlang::syms(c(...))
  for (s in cands) if (rlang::as_string(s) %in% names(df)) return(rlang::as_string(s))
  stop("Could not find the corresponding column: ", paste(c(...), collapse = ", "))
}

chr_col   <- std_col(qtl_summary, "CHROM", "chr", "chromosome")
start_col <- std_col(qtl_summary, "startPos", "start.pos", "start", "startBP", "startPOS")
end_col   <- std_col(qtl_summary, "endPos",   "end.pos",   "end",   "endBP",   "endPOS")

# Smoothing by number of SNPs - running mean over k points
smooth_by_snp <- function(x, k = 100) {
  zoo::rollapply(x, width = k, FUN = mean, align = "center", fill = NA)
}

# Function that plots two vertically stacked "tracks" for one QTL:
# top = |deltaSNP| (smoothed), bottom = SNP indexes of both bulks (smoothed, colored)
plot_qtl_tracks <- function(dat, chr, start_pos, end_pos,
                            low_name  = lowBulk,
                            high_name = highBulk,
                            smooth_k = 100,
                            threshold = 0.04,
                            return_intervals = FALSE,
                            return_peak = FALSE) {

  subdat <- dat %>%
    dplyr::filter(CHROM == chr, POS >= start_pos, POS <= end_pos) %>%
    dplyr::arrange(POS)

  subdat <- subdat %>%
    mutate(
      deltaSNP_abs             = abs(deltaSNP),
      deltaSNP_smooth          = smooth_by_snp(deltaSNP_abs, k = smooth_k),
      SNPindex_lowbulk_smooth  = smooth_by_snp(SNPindex_lowbulk, k = smooth_k),
      SNPindex_highbulk_smooth = smooth_by_snp(SNPindex_highbulk, k = smooth_k)
    )

  # intervals above threshold (computed from |deltaSNP|)
  runs <- rle(subdat$deltaSNP_smooth > threshold)
  ends <- cumsum(runs$lengths)
  starts <- c(1, head(ends, -1) + 1)
  intervals <- data.frame(
    CHROM    = chr,
    startPOS = subdat$POS[starts],
    endPOS   = subdat$POS[ends],
    above    = runs$values
  )
  intervals_above <- intervals %>% filter(above)

  # position of the maximum (highest peak) in absolute deltaSNP
  peak_idx <- which.max(subdat$deltaSNP_smooth)
  peak_pos <- subdat$POS[peak_idx]
  peak_val <- subdat$deltaSNP_smooth[peak_idx]
  peak_df <- data.frame(CHROM = chr, peakPOS = peak_pos, peakValue = peak_val)

  if (return_peak) return(peak_df)
  if (return_intervals) return(intervals_above)

  # Panel 1: |deltaSNP|
  p_delta <- ggplot(subdat, aes(x = POS, y = deltaSNP_smooth)) +
    geom_line(size = 0.4) +
    labs(title = paste0(chr, ": ", scales::comma(start_pos), "-", scales::comma(end_pos),
                        " (|deltaSNP|, smoothed ", smooth_k, " SNP)"),
         x = NULL,
         y = expression("|" * Delta * "SNP|")) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          axis.title.x = element_blank())

  # Panel 2: SNP index of both bulks
  p_snpidx <- ggplot(subdat, aes(x = POS)) +
    geom_line(aes(y = SNPindex_lowbulk_smooth,  color = low_name),  size = 0.5) +
    geom_line(aes(y = SNPindex_highbulk_smooth, color = high_name), size = 0.5) +
    scale_color_manual(values = c("steelblue", "darkorange"),
                       breaks = c(low_name, high_name),
                       name = NULL) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(title = paste0("SNP index (smoothed ", smooth_k, " SNP)"),
         x = "Genomic position",
         y = "SNP index") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  cowplot::plot_grid(p_delta, p_snpidx, ncol = 1, align = "v", rel_heights = c(1,1))
}

# --- main part ---
all_intervals <- list()
all_peaks <- list()

pdf(paste0("QTL_tracks_smoothed_", tag, "_ABS.pdf"), width = 10, height = 6)

for (i in seq_len(nrow(qtl_summary))) {
  chr_i   <- qtl_summary[[chr_col]][i]
  start_i <- as.numeric(qtl_summary[[start_col]][i])
  end_i   <- as.numeric(qtl_summary[[end_col]][i])

  if (!is.na(start_i) && !is.na(end_i) && start_i > end_i) {
    tmp <- start_i; start_i <- end_i; end_i <- tmp
  }

  # intervals
  intervals_i <- plot_qtl_tracks(df_filt, chr_i, start_i, end_i,
                                 low_name = lowBulk, high_name = highBulk,
                                 return_intervals = TRUE)
  if (nrow(intervals_i) > 0) all_intervals[[length(all_intervals)+1]] <- intervals_i

  # peak position (not plotted)
  peak_i <- plot_qtl_tracks(df_filt, chr_i, start_i, end_i,
                            low_name = lowBulk, high_name = highBulk,
                            return_peak = TRUE)
  all_peaks[[length(all_peaks)+1]] <- peak_i

  # plot
  print(
    plot_qtl_tracks(df_filt, chr_i, start_i, end_i,
                    low_name = lowBulk, high_name = highBulk,
                    smooth_k = 100, threshold = 0.04,
                    return_intervals = FALSE)
  )
}

dev.off()

# Save intervals above threshold and peak maxima
if (length(all_intervals) > 0) {
  all_intervals_df <- dplyr::bind_rows(all_intervals)
  readr::write_csv(all_intervals_df, paste0("QTL_deltaSNP_intervals_above0.04_", tag, "_ABS.csv"))
}

if (length(all_peaks) > 0) {
  all_peaks_df <- dplyr::bind_rows(all_peaks)
  readr::write_csv(all_peaks_df, paste0("QTL_deltaSNP_peak_positions_", tag, "_ABS.csv"))
  cat("\nTop |deltaSNP| peaks per QTL region:\n")
  print(all_peaks_df)
}

##### Gene Annotation and Candidate Matching #####

library(ape)

# Input files (see the header of this script and the README)
gtf_path   <- "input/Arabidopsis_arenosa_genome.annotation.gtf.gz"
ortho_path <- "input/Arabidopsis_arenosa_At_orthologues_and_GO.tsv"

# FST candidate genes: the divergence scan assigns each window to a gene using the gene
# body OR an upstream/downstream regulatory region up to 2 kb away (clipped to the
# neighbouring gene if that is closer). To compare the two sets on the same footing, the
# same gene + regulatory-region extension is built below over the FULL gene annotation and
# used as the SNP-to-gene overlap table, so that a significant SNP in a gene's flank counts
# as a hit for that gene, as it does on the FST side.
fst_candidates_path <- "input/fst_candidates.tsv"
PCT_THRESHOLD <- 99  # top 1% (fstwc.percentile >= 99)

# A. arenosa gene annotation (AUGUSTUS GTF). The attribute column contains just the bare
# gene ID (e.g. "jg10308"). Only the feature type "gene" is used.
annotation_raw <- read.gff(gzfile(gtf_path))
genes_gff <- annotation_raw[annotation_raw$type == "gene", ]
genes_gff$gene <- genes_gff$attributes
genes_gff <- as.data.table(genes_gff)
colnames(genes_gff)[1] <- "scaff"
genes_gff <- genes_gff[scaff %in% Chroms]
# Sort by scaffold/start/end - "previous"/"next" gene below means row-adjacent
# within the same scaffold under this order.
setkey(genes_gff, scaff, start, end)

## Upstream region: 2kb before the gene start, or up to the end of the
## previous gene if that's closer (chr_previous guard prevents crossing a
## scaffold boundary). Labelled "upstream_region"/"downstream_region"
## depending on strand - a region lying before start_of_this is biologically
## downstream for a minus-strand gene.
upstr_gff <- copy(genes_gff)
upstr_gff$start <- mapply(
  FUN = function(chr_this, chr_previous, start_of_this, stop_of_previous) {
    if (chr_this == chr_previous & start_of_this > stop_of_previous) {
      max(start_of_this - 2000, stop_of_previous + 1)
    } else {
      max(start_of_this - 2000, 0 + 1)
    }
  },
  chr_this         = genes_gff$scaff,
  chr_previous     = c("none", head(genes_gff$scaff, -1)),
  start_of_this    = genes_gff$start,
  stop_of_previous = c(0, head(genes_gff$end, -1))
)
upstr_gff$end  <- genes_gff$start - 1
upstr_gff$type <- ifelse(genes_gff$strand == "+", "upstream_region", "downstream_region")

## Downstream region: 2kb after the gene end, or up to the start of the next
## gene if that's closer. Symmetric counterpart of the upstream block above.
downstr_gff <- copy(genes_gff)
downstr_gff$end <- mapply(
  FUN = function(chr_this, chr_next, end_of_this, start_of_next) {
    if (chr_this == chr_next & end_of_this < start_of_next) {
      min(end_of_this + 2000, start_of_next - 1)
    } else {
      end_of_this + 2000
    }
  },
  chr_this      = genes_gff$scaff,
  chr_next      = c(tail(genes_gff$scaff, -1), "none"),
  end_of_this   = genes_gff$end,
  start_of_next = c(tail(genes_gff$start, -1), 0)
)
downstr_gff$start <- genes_gff$end + 1
downstr_gff$type  <- ifelse(genes_gff$strand == "+", "downstream_region", "upstream_region")

# Combine gene body + its regulatory regions into one SNP-to-gene overlap
# table. A SNP can now land in a gene's body OR its regulatory flank and
# still count as hitting that gene - matching how fst_candidates_path assigns
# its windows to genes.
annotation2 <- rbind(genes_gff, upstr_gff, downstr_gff, fill = TRUE)
annotation2 <- annotation2[end != 0]  # drop 0-length regions (gene starting at scaffold position 1)
setkey(annotation2, scaff, start, end)

# Gene-body-only table, kept separately for gene_coords below (which needs
# the gene's own coordinates, not a regulatory region's - annotation2 now has
# up to 3 rows per gene: body, upstream, downstream).
annotation_genes <- genes_gff

all_sig_snps <- read.csv(paste0("BSA_QTLseqr_significant_SNPs_", tag, "_FDR0.05.csv"), header = TRUE)

######## SNP-level gene annotation - take all significant SNPs, assign them to genes
######## (gene body or regulatory flank, see annotation2 above), count how many
######## significant SNPs fall in/near each gene and pull the gene's stats
all_sig_snps <- all_sig_snps %>%
  rename(scaff = CHROM, start = POS) %>%
  mutate(end = start)

setDT(all_sig_snps)
all_sig_snps <- foverlaps(all_sig_snps, annotation2, type = "any")
all_sig_snps <- subset(all_sig_snps, !is.na(gene))

# assign gene-level stats from the SNPs that overlap it
snps_to_gene <- all_sig_snps %>%
  filter(!is.na(gene)) %>%
  # guard against the same SNP overlapping the same gene more than once (e.g.
  # matching both the gene body row and, in principle, a regulatory-region row)
  distinct(gene, scaff, i.start, .keep_all = TRUE) %>%
  group_by(gene) %>%
  summarise(
    pocet_snps = n(),

    deltaSNP_max      = max(deltaSNP, na.rm = TRUE),
    deltaSNP_mean     = mean(deltaSNP, na.rm = TRUE),
    deltaSNP_abs_max  = max(abs(deltaSNP), na.rm = TRUE),
    deltaSNP_abs_mean = mean(abs(deltaSNP), na.rm = TRUE),

    G_max          = max(G, na.rm = TRUE),
    G_mean         = mean(G, na.rm = TRUE),
    Gprime_max     = max(Gprime, na.rm = TRUE),
    Gprime_mean    = mean(Gprime, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  arrange(desc(Gprime_max))

# A. arenosa gene -> A. thaliana ortholog (the ortholog table is keyed by the same gene IDs
# as the GTF). Matching of the QTL and FST gene sets is done on A. arenosa gene IDs;
# the ortholog columns (single_Arabidopsis_thaliana, UniProt_Protein.names) are used only
# for gene descriptions in the candidate table and plots.
orthologs <- read.delim(ortho_path, header = TRUE, stringsAsFactors = FALSE)
snps_to_gene <- snps_to_gene %>%
  left_join(orthologs, by = c("gene" = "Arabidopsis_arenosa"))

# how many candidate genes have a single A. thaliana ortholog (for reference only)
sum(!is.na(snps_to_gene$single_Arabidopsis_thaliana) & snps_to_gene$single_Arabidopsis_thaliana != "")

# attach which QTL each gene falls in, and its genomic coordinates. Coordinates come from
# annotation_genes (gene body only), not annotation2, which also contains the
# regulatory-region rows.
gene_coords <- annotation_genes %>%
  as.data.frame() %>%
  select(gene, scaff, gene_start = start, gene_end = end) %>%
  distinct(gene, .keep_all = TRUE)

gene_qtl <- all_sig_snps %>%
  as.data.frame() %>%
  filter(!is.na(gene)) %>%
  group_by(gene) %>%
  summarise(qtl = unique(qtl)[1], .groups = "drop")

snps_to_gene <- snps_to_gene %>%
  left_join(gene_qtl, by = "gene") %>%
  left_join(gene_coords, by = "gene") %>%
  relocate(qtl, scaff, gene_start, gene_end, .after = pocet_snps)

write.csv(snps_to_gene, paste0("candidate_genes_SNPapproach_BSA_", tag, ".csv"), row.names = FALSE)

# Compare the QTL genes with the FST candidate genes of the BOR-VLA parental population
# pair (design as in https://doi.org/10.1038/s41467-021-25256-5), matched by A. arenosa
# gene ID. The "BOR_VLA_intersection" column and the "BOR-VLA" file names refer to this pair.
FST_candidates <- read.delim(fst_candidates_path, header = TRUE, stringsAsFactors = FALSE)
FST_candidates <- FST_candidates[FST_candidates$fstwc.percentile >= PCT_THRESHOLD, ]
# Restrict to the same 8 scaffolds the QTL analysis is run on (Chroms, defined
# above) - the candidate file also carries a handful of windows on small
# unplaced scaffolds (scaffold_9 and up) that snps_to_gene can never contain.
FST_candidates <- FST_candidates[FST_candidates$chr %in% Chroms, ]
FST_candidate_genes <- unique(FST_candidates$gene_ID)
FST_candidate_genes <- FST_candidate_genes[!is.na(FST_candidate_genes)]

snps_to_gene$BOR_VLA_intersection <- ifelse(snps_to_gene$gene %in% FST_candidate_genes, 1, 0)
matched_candidates_single_SNPs <- snps_to_gene[snps_to_gene$BOR_VLA_intersection == 1, ]
write.csv(matched_candidates_single_SNPs, paste0("matched_candidates_SNPapproach_BSA_", tag, ".csv"), row.names = FALSE)

##### GO enrichment: run in a separate script (see the GO enrichment folder) #####

#################### Plots for matched candidate genes ####################
# close-ups of |deltaSNP| and SNPindex around matched candidates (auto-scaled range + peak table)

# === Settings ===
flank_size <- 30000
smooth_k <- 100

# === Safer smoothing ===
smooth_by_snp <- function(x, k = 100) {
  if (length(x) < k) {
    return(rep(mean(x, na.rm = TRUE), length(x)))
  }
  zoo::rollapply(x, width = k, FUN = mean, align = "center", fill = NA)
}

# === Clip interval to what data is actually available ===
clip_interval <- function(pos, start_pos, end_pos, flank_bp) {
  min_pos <- min(pos, na.rm = TRUE)
  max_pos <- max(pos, na.rm = TRUE)
  start_clipped <- max(start_pos - flank_bp, min_pos)
  end_clipped   <- min(end_pos + flank_bp,   max_pos)
  return(c(start_clipped, end_clipped))
}

# === Plotting + computation ===
plot_gene_tracks <- function(dat, chr, start_pos, end_pos,
                             flank_bp = flank_size,
                             low_name  = lowBulk,
                             high_name = highBulk,
                             smooth_k = smooth_k,
                             gene_id = NULL,
                             short_desc = NULL) {

  pos_chr <- dat$POS[dat$CHROM == chr]
  if (length(pos_chr) == 0) return(NULL)
  bounds <- clip_interval(pos_chr, start_pos, end_pos, flank_bp)

  sub <- dat %>%
    filter(CHROM == chr, POS >= bounds[1], POS <= bounds[2]) %>%
    arrange(POS)

  if (nrow(sub) == 0) return(NULL)

  sub <- sub %>%
    mutate(
      deltaSNP_abs             = abs(deltaSNP),
      deltaSNP_smooth          = smooth_by_snp(deltaSNP_abs,      k = smooth_k),
      SNPindex_lowbulk_smooth  = smooth_by_snp(SNPindex_lowbulk,  k = smooth_k),
      SNPindex_highbulk_smooth = smooth_by_snp(SNPindex_highbulk, k = smooth_k)
    )

  valid_idx <- which(!is.na(sub$deltaSNP_smooth))
  shown_start <- min(sub$POS[valid_idx])
  shown_end   <- max(sub$POS[valid_idx])

  shown_left  <- round((start_pos - shown_start) / 1000, 1)
  shown_right <- round((shown_end - end_pos) / 1000, 1)
  range_label <- paste0("(shown: -", shown_left, " kb / +", shown_right, " kb)")

  peak_idx <- which.max(sub$deltaSNP_smooth)
  peak_pos <- sub$POS[peak_idx]
  peak_val <- sub$deltaSNP_smooth[peak_idx]

  # Panel 1: |deltaSNP|
  p_delta <- ggplot(sub, aes(x = POS, y = deltaSNP_smooth)) +
    geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = -Inf, ymax = Inf),
              fill = "grey50", alpha = 0.6, inherit.aes = FALSE) +
    geom_line(size = 0.4, color = "black", na.rm = TRUE) +
    coord_cartesian(xlim = c(shown_start, shown_end)) +
    labs(title = paste0(gene_id, " | ", chr, ": ",
                        scales::comma(start_pos), "-", scales::comma(end_pos),
                        " ", range_label),
         subtitle = short_desc,
         x = NULL, y = expression("|" * Delta * "SNP|")) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9, face = "italic"))

  # Panel 2: SNP index
  p_snpidx <- ggplot(sub, aes(x = POS)) +
    geom_rect(aes(xmin = start_pos, xmax = end_pos, ymin = -Inf, ymax = Inf),
              fill = "grey50", alpha = 0.6, inherit.aes = FALSE) +
    geom_line(aes(y = SNPindex_lowbulk_smooth,  color = low_name),
              size = 0.5, na.rm = TRUE) +
    geom_line(aes(y = SNPindex_highbulk_smooth, color = high_name),
              size = 0.5, na.rm = TRUE) +
    coord_cartesian(ylim = c(0, 1), xlim = c(shown_start, shown_end)) +
    scale_color_manual(values = c("steelblue", "darkorange"),
                       breaks = c(low_name, high_name), name = NULL) +
    labs(x = "Genomic position", y = "SNP index") +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"),
          legend.position = "top")

  plot_obj <- cowplot::plot_grid(p_delta, p_snpidx, ncol = 1, align = "v", rel_heights = c(1,1))

  data.frame(
    gene_id = gene_id,
    CHROM = chr,
    start = start_pos,
    end = end_pos,
    shown_start = shown_start,
    shown_end = shown_end,
    shown_left_kb = shown_left,
    shown_right_kb = shown_right,
    peakPOS = peak_pos,
    peakValue = round(peak_val, 5)
  ) -> peak_info

  return(list(plot = plot_obj, info = peak_info))
}

# === Export PDF and CSV ===
# NOTE: this loop runs over matched_candidates_single_SNPs (the SNP-level
# candidate/FST-overlap table actually computed above). Gene descriptions come
# from the ortholog table's UniProt_Protein.names column.
pdf(paste0("Matched_candidate_genes_tracks_", tag, "_ABS_30kbflanking.pdf"), width = 10, height = 6)
peak_summary <- list()

for (i in seq_len(nrow(matched_candidates_single_SNPs))) {
  chr_i   <- matched_candidates_single_SNPs$scaff[i]
  start_i <- as.numeric(matched_candidates_single_SNPs$gene_start[i])
  end_i   <- as.numeric(matched_candidates_single_SNPs$gene_end[i])
  gene_id <- matched_candidates_single_SNPs$gene[i]
  desc_i  <- matched_candidates_single_SNPs$UniProt_Protein.names[i]

  if (!is.na(start_i) && !is.na(end_i)) {
    if (start_i > end_i) { tmp <- start_i; start_i <- end_i; end_i <- tmp }

    message("Plotting ", gene_id, " (", chr_i, ":", start_i, "-", end_i, ")")

    result <- plot_gene_tracks(df_filt,
                               chr = chr_i, start_pos = start_i, end_pos = end_i,
                               flank_bp = flank_size,
                               low_name = lowBulk, high_name = highBulk,
                               smooth_k = smooth_k,
                               gene_id = gene_id,
                               short_desc = desc_i)

    if (!is.null(result)) {
      print(result$plot)
      peak_summary[[length(peak_summary) + 1]] <- result$info
    }
  }
}

dev.off()

if (length(peak_summary) > 0) {
  peak_summary_df <- bind_rows(peak_summary)
  write_csv(peak_summary_df, paste0("Matched_candidate_genes_tracks_", tag, "_ABS_30kbflanking.csv"))
  message("Table saved as 'Matched_candidate_genes_tracks_", tag, "_ABS_30kbflanking.csv'")
}

############### Superexact test ###############

# make a vector of gene names for our QTL candidate genes and the FST candidate
# genes, without duplicates. Both sets are arenosa gene IDs (matching how the
# intersection above is computed), not A. thaliana ortholog codes.
clean_set <- function(x){
  x <- trimws(x)
  x <- x[!is.na(x) & x != ""]
  unique(x)
}

BSA_QTL_final_set <- clean_set(snps_to_gene$gene)
FST_final_set      <- clean_set(FST_candidate_genes)
length(BSA_QTL_final_set)
length(FST_final_set)

# For the SuperExact test, build a list of the sets to compare
all_res <- list(BSA = BSA_QTL_final_set,
                 FST = FST_final_set)

# Background population size (n) for the test: the total number of annotated A. arenosa
# genes on the same 8 scaffolds the QTL analysis is run on (annotation_genes, gene body
# only, already restricted to Chroms).
total <- length(unique(annotation_genes$gene))

library("SuperExactTest")
res <- supertest(all_res, n = total)
summary(res)
write.csv(summary(res)$Table, file = paste0("Fisher_exact_test_BSAQTLgenes_BOR-VLA_", tag, ".csv"), row.names = FALSE)

pdf(paste0("Superexact_test_BSAQTLgenes_BOR-VLA_", tag, ".pdf"), width = 10, height = 7)
plot(res, Layout="landscape", degree=1:2, sort.by="size", margin=c(0.5,6,2,3), keep.empty.intersections=F,
     cex=0.7, cex.lab= 0.7, legend.col=2, legend.text.cex=0.5, color.scale.cex=0.5,color.on='black', color.off='#EEEEEE',overlap.size.cex=0.5, track.area.range=0.4, bar.area.range=0.2)
dev.off()

#################################################################################

#Venn diagramme
library(ggvenn)

venn_list <- list(
  "BSA" = BSA_QTL_final_set,
  "FST" = FST_final_set
)

pdf(paste0("Venn_BSA_FST_", tag, ".pdf"), width = 6, height = 6)
ggvenn(
  venn_list,
  fill_color = c("#ffc300", "#003566"),
  stroke_size = 0.5,
  set_name_size = 4,
  text_size = 5,
  show_percentage = FALSE
)
dev.off()
