# Step 3: F_ST outlier windows -> candidate genes.
#
# Reads the 1-kb window scan from 02_piawka_fst_scan.sh, keeps windows with
# >= MIN_LINES SNPs and <= MAX_MISS missing genotypes in both populations,
# ranks windows by Weir & Cockerham F_ST (fstwc), and assigns the windows to genes.
# A window is assigned to a gene when it overlaps the gene body or its flanks:
# 2 kb upstream of the start and 2 kb downstream of the end (measured from the start
# and stop codon because the annotation has no UTRs), shortened where the flank
# would reach the neighbouring gene.
#
# Output:
#   fst_candidates.tsv         one row per window x gene region (gene body / upstream /
#                              downstream) for windows above TOP_QUANTILE_OUT (default 0.95);
#                              script 03_qtl_analysis keeps fstwc.percentile >= 99 (top 1 %)
#   fst_candidate_genes_top1pct.tsv   unique genes for the top 1 % windows
#
# Usage: Rscript 03_fst_candidates.R Fst_scan_VLA_BOR.bed.gz annotation.gtf output_dir
# Requires: data.table, R.utils (for reading .gz)

suppressPackageStartupMessages(library(data.table))

LINES_THRESHOLD  <- 10    # minimum number of SNPs per window in both populations
MISS_THRESHOLD   <- 0.5   # maximum share of missing genotypes in both populations
FLANK            <- 2000  # bp
TOP_QUANTILE_OUT <- 0.95  # windows above this quantile are written to fst_candidates.tsv
TOP_PERCENTILE   <- 99    # top 1 % of windows
MAIN_SCAFFOLDS   <- paste0("scaffold_", 1:8)   # the eight assembled chromosomes

## --- scan -> one row per window with all statistics ---------------------------------
read_scan <- function(scan_file) {
  sc <- fread(scan_file, header = TRUE)
  setnames(sc, sub("#", "", names(sc), fixed = TRUE))

  ## window skeleton and between-population statistics
  wide <- unique(sc[, 1:4])
  wide[, wCenter := as.integer(rowMeans(.SD)), .SDcols = c("start", "end")]
  wide[, chrNum := as.integer(sub("scaffold_", "", chr))]
  for (v in c("fst", "fstwc", "afd")) {
    cols <- if (v == "fst") c(1:6, 8) else c(1:4, 8)
    one <- sc[stat == v, ..cols]
    wide <- merge(wide, one, by = c("chr", "start", "end", "locus"), all.x = TRUE)
    setnames(wide, "value", v)
  }
  pop1 <- wide$pop1[!is.na(wide$pop1)][1]
  pop2 <- wide$pop2[!is.na(wide$pop1)][1]

  ## within-population statistics
  for (v in c("lines", "miss", "segr")) {
    for (k in 1:2) {
      pp <- c(pop1, pop2)[k]
      one <- sc[stat == v & pop1 == pp, c(1:4, 8)]
      wide <- merge(wide, one, by = c("chr", "start", "end", "locus"), all.x = TRUE)
      setnames(wide, "value", paste0(v, ".pop", k))
    }
  }
  wide
}

## --- gene bodies + upstream / downstream regions --------------------------------------
build_gene_regions <- function(gtf_file, flank = FLANK) {
  ## the BRAKER annotation contains many "#" comment lines, which read.table skips
  gff <- as.data.table(read.table(gtf_file, header = FALSE, sep = "\t", quote = "",
                                  comment.char = "#", stringsAsFactors = FALSE))
  genes <- gff[grepl("gene", V3)]                 # genes and pseudogenes
  genes[, gene_ID := V9]
  setkey(genes, V1, V4, V5)
  n <- nrow(genes)

  ## upstream region: up to `flank` bp before the start, or up to the end of the previous gene
  same_prev <- c(FALSE, genes$V1[-1] == genes$V1[-n])
  prev_end  <- c(0L, genes$V5[-n])
  up <- copy(genes)
  up[, V4 := ifelse(same_prev & genes$V4 > prev_end,
                    pmax(genes$V4 - flank, prev_end + 1L),
                    pmax(genes$V4 - flank, 1L))]
  up[, V5 := genes$V4 - 1L]
  up[, V3 := ifelse(genes$V7 == "+", "upstream_region", "downstream_region")]

  ## downstream region: up to `flank` bp after the end, or up to the start of the next gene
  same_next  <- c(genes$V1[-n] == genes$V1[-1], FALSE)
  next_start <- c(genes$V4[-1], 0L)
  down <- copy(genes)
  down[, V5 := ifelse(same_next & genes$V5 < next_start,
                      pmin(genes$V5 + flank, next_start - 1L),
                      genes$V5 + flank)]
  down[, V4 := genes$V5 + 1L]
  down[, V3 := ifelse(genes$V7 == "+", "downstream_region", "upstream_region")]

  regions <- rbind(genes, up, down)
  regions <- regions[V5 != 0]                      # gene starting at position 1 has no upstream region
  setkey(regions, V1, V4, V5)
  regions
}

## --- windows -> genes ------------------------------------------------------------------
annotate_windows <- function(windows, regions) {
  setkey(windows, chr, start, end)
  foverlaps(x = windows, y = regions, by.x = c("chr", "start", "end"),
            by.y = c("V1", "V4", "V5"), type = "any")
}

OUT_COLUMNS <- c("gene_ID", "chr", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9",
                 "start", "end", "locus", "wCenter", "chrNum", "pop1", "pop2",
                 "fst", "fstwc", "afd", "lines.pop1", "lines.pop2", "miss.pop1", "miss.pop2",
                 "segr.pop1", "segr.pop2", "fstwc.percentile")

run_scan_to_candidates <- function(scan_file, gtf_file, out_dir) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  win <- read_scan(scan_file)

  ## windows with enough SNPs and low missingness in both populations
  win <- win[lines.pop1 >= LINES_THRESHOLD & lines.pop2 >= LINES_THRESHOLD &
             miss.pop1 <= MISS_THRESHOLD & miss.pop2 <= MISS_THRESHOLD]
  win[, fstwc.percentile := rank(fstwc, ties.method = "min", na.last = FALSE) / .N * 100]
  cat("windows after filtering:", nrow(win), "\n")

  regions <- build_gene_regions(gtf_file)
  all_windows <- annotate_windows(win, regions)

  top <- all_windows[fstwc.percentile >= TOP_QUANTILE_OUT * 100]
  top <- top[order(-fstwc)]
  fwrite(top[, ..OUT_COLUMNS], file.path(out_dir, "fst_candidates.tsv"), sep = "\t", quote = FALSE, na = "NA")

  top1 <- top[fstwc.percentile >= TOP_PERCENTILE & chr %in% MAIN_SCAFFOLDS & !is.na(gene_ID)]
  genes <- sort(unique(top1$gene_ID))
  fwrite(data.table(gene_ID = genes), file.path(out_dir, "fst_candidate_genes_top1pct.tsv"), sep = "\t")
  cat("top 1 % windows (unique):", uniqueN(all_windows[fstwc.percentile >= TOP_PERCENTILE, .(chr, start, end)]),
      "\ncandidate genes on scaffolds 1-8:", length(genes), "\n")
  invisible(top)
}

if (!interactive() && sys.nframe() == 0) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) != 3) stop("Usage: Rscript 03_fst_candidates.R Fst_scan.bed.gz annotation.gtf output_dir")
  run_scan_to_candidates(args[1], args[2], args[3])
}
