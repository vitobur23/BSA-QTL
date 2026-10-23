####################################################################################
# GO enrichment of the genes shared by the BSA QTL regions and the FST divergence-scan
# candidates (A. arenosa)
#
# Input: the 61 genes in the overlap of the BSA QTL genes and the FST candidates
# (output of ../03_qtl_analysis). Enrichment is tested with topGO (Fisher's exact test,
# "elim" algorithm) against all annotated A. arenosa genes on the eight assembled
# scaffolds.
#
# Run from this folder (setwd() below, or in RStudio: Session > Set Working Directory >
# To Source File Location). All paths are relative to this folder and all output files
# are written here.
#
# GO source: GO terms of each A. arenosa gene were transferred from its A. thaliana
# ortholog when the ortholog table was built
# (https://github.com/mduchoslav/Brassicaceae_orthology). The enrichment itself runs
# entirely at the A. arenosa gene level. The table is keyed by the same jgXXXXX gene IDs
# as the gene annotation and the candidate table.
#
# The test is run ONLY on the overlap set, not on all genes inside the QTL intervals.
# The QTL intervals are wide and most genes in them are there through linkage, so an
# enrichment test on all of them would describe the intervals rather than the trait.
#
# Outputs: significant terms per ontology (BP, MF, CC), merged BP+MF table and plots
# (dotplot: colour = p-value on a log scale, point size = gene count; bar plot of fold
# enrichment), and a table of the candidate genes annotated to each significant term.
####################################################################################

# setwd("/path/to/04_go_enrichment")


########
#0. Packages
##########

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("topGO", quietly = TRUE)) {
  BiocManager::install("topGO", ask = FALSE, update = FALSE)
  # GO.db can fail on the first pass if its own dependencies are still being
  # compiled in the same call - installing it again afterwards fixes that.
  if (!requireNamespace("GO.db", quietly = TRUE)) {
    BiocManager::install("GO.db", ask = FALSE, update = FALSE)
  }
}
for (p in c("ggplot2", "dplyr", "stringr", "viridisLite")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}

library(topGO)
library(ggplot2)
library(dplyr)
library(stringr)


########
#0b. Settings and input paths
##########

# tag carried by every output file
tag   <- "overlap61_ws1000000_DP4_arenosa_top1pct"
title <- "BSA QTL x FST overlap genes (n = 61)"

# A. arenosa gene -> A. thaliana ortholog + GO annotation
# (the same file the QTL analysis script uses for gene descriptions)
ortho_path <- "../03_qtl_analysis/input/Arabidopsis_arenosa_At_orthologues_and_GO.tsv"

# Gene annotation (AUGUSTUS GTF) - only used to build the gene universe, i.e. which genes
# sit on the 8 assembled scaffolds the QTL analysis was run on
gtf_path <- "../03_qtl_analysis/input/Arabidopsis_arenosa_genome.annotation.gtf.gz"

# The 61 overlap genes (BSA QTL genes x FST / BOR-VLA candidates), output of the QTL
# analysis script (Supplementary Data S4)
cand_path <- "../03_qtl_analysis/matched_candidates_SNPapproach_BSA_ws1000000_DP4_arenosa_top1pct.csv"

# Gene universe (= background) for the enrichment test.
#   TRUE  - only genes on scaffold_1..8, i.e. exactly the genes that could have come out
#           as candidates. This matches the background used for the Fisher / SuperExact
#           test in the QTL analysis script (33,275 genes).
#   FALSE - every gene in the ortholog table.
RESTRICT_UNIVERSE_TO_QTL_SCAFFOLDS <- TRUE
Chroms <- paste0("scaffold_", 1:8)

# Upper limit on how many genes a GO term may be annotated to. Broad terms (annotated to
# more than 500 genes) are excluded from the merged BP+MF table and plots, as are terms
# represented by a single candidate gene. Set to Inf to keep broad terms.
MAX_ANNOTATED <- 500


########
#1. Preparation of the GO mapping library
##########

# quote = "", comment.char = "" is needed to avoid mess resulting from some gene names
ortho_df <- read.table(file = ortho_path, header = TRUE, sep = "\t",
                       quote = "", comment.char = "", na.strings = "")
dim(ortho_df)
colnames(ortho_df)
# Arabidopsis_arenosa | prot_length | Arabidopsis_thaliana |
# single_Arabidopsis_thaliana | homologue_type | GO_term_IDs | UniProt_Protein.names

# genes with GO terms
sum(!is.na(ortho_df$GO_term_IDs))

# write just the GO mappings: gene ID (col 1) + GO terms (col 6)
write.table(x = ortho_df[, c(1, 6)],
            file = "GO_mappings_Arabidopsis_arenosa.tsv",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# use topGO to read the GO mappings
geneID2GO <- readMappings("GO_mappings_Arabidopsis_arenosa.tsv", sep = "\t", IDsep = ";")

# define the gene universe
geneUniverse <- names(geneID2GO)
length(geneUniverse)

if (RESTRICT_UNIVERSE_TO_QTL_SCAFFOLDS) {
  # gene rows of the AUGUSTUS GTF; the attribute column (V9) is the bare gene ID
  gtf <- read.table(gzfile(gtf_path), sep = "\t", header = FALSE, quote = "",
                    comment.char = "#", fill = TRUE, stringsAsFactors = FALSE)
  gtf <- gtf[gtf$V3 == "gene" & gtf$V1 %in% Chroms, ]
  genes_on_chroms <- unique(gtf$V9)
  cat("genes on the", length(Chroms), "assembled scaffolds:", length(genes_on_chroms), "\n")

  geneID2GO    <- geneID2GO[names(geneID2GO) %in% genes_on_chroms]
  geneUniverse <- names(geneID2GO)
}
cat("gene universe size:", length(geneUniverse), "\n")


########
#1b. Candidate genes
##########

cand <- read.table(cand_path, header = TRUE, sep = ",", quote = "\"",
                   comment.char = "", stringsAsFactors = FALSE)
nrow(cand)                        # 61 overlap genes
cand_genes <- unique(cand$gene)
length(cand_genes)

# sanity check - every candidate should be findable in the universe
cat("candidates not in the gene universe:", sum(!cand_genes %in% geneUniverse), "\n")

# presence/absence vector for the whole gene universe
presence.cand <- factor(as.integer(geneUniverse %in% cand_genes))
names(presence.cand) <- geneUniverse

# check the number of genes present
sum(presence.cand == 1)

int.genes <- presence.cand


# Re-definition of a topGO method to change the rounding of "Expected" (5 digits
# instead of 3), so that small expected values are not rounded to zero
setMethod("termStat",
					signature(object = "topGOdata", whichGO = "character"),
					function(object, whichGO) {

						if(length(whichGO) == 0) {
							return(data.frame(Annotated = integer(),
																Significant = integer(),
																Expected = numeric()))
						}

						x <- .genesInNode(graph(object), whichGO)

						anno <- sapply(x, length)

						sGenes <- sigGenes(object)
						sig <- sapply(x, function(e) length(intersect(e, sGenes)))

						expect <- numSigGenes(object) / length(genes(object))
						expect <- round(anno * expect, 5) # rounding changed from 3 to 5 digits

						return(data.frame(Annotated = anno,
															Significant = sig,
															Expected = expect,
															row.names = names(x)))
					})


########
#2. Helpers
##########

# GenTable returns p-values as character ("< 1e-30" for very small ones), so always
# go through this before comparing or plotting them.
p_num <- function(x) suppressWarnings(as.numeric(sub("^\\s*<\\s*", "", x)))

# fold enrichment (FE) = Significant / Expected
add_FE <- function(x) {
  x$Annotated   <- as.numeric(x$Annotated)
  x$Significant <- as.numeric(x$Significant)
  x$Expected    <- as.numeric(x$Expected)
  x$Expected[x$Expected == 0] <- 0.001   # adjustment as you cannot divide by 0
  x$FE <- x$Significant / x$Expected
  x$elimFisher_num <- p_num(x$elimFisher)
  x
}

run_go <- function(ontology, order_by = "elimFisher") {
  go.obj <- new("topGOdata", ontology = ontology, allGenes = int.genes,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  # enrichment tests: elim (more conservative, accounts for the GO hierarchy) and classic
  resultsFe <- runTest(go.obj, algorithm = "elim",    statistic = "fisher")
  resultsFc <- runTest(go.obj, algorithm = "classic", statistic = "fisher")
  # topNodes must not exceed the number of GO terms actually in the DAG - the CC
  # ontology has fewer than 1000 here, and asking for more makes GenTable fail
  # with "arguments imply differing number of rows".
  n_top <- min(1000, length(usedGO(go.obj)))
  # numChar = 1000 keeps the FULL term name. GenTable's default is 40, which cuts
  # names off with "...", which makes the plot labels unreadable.
  tab <- GenTable(go.obj, classicFisher = resultsFc, elimFisher = resultsFe,
                  orderBy = order_by, ranksOf = order_by, topNodes = n_top,
                  numChar = 1000)
  list(obj = go.obj, tab = tab)
}

# Size-legend breaks for the dotplot. The gene counts here are small integers
# (typically 1-4), so fixed breaks such as c(1,3,5,10,20) would leave a single bubble
# in the legend. Take up to 3 integer breaks spanning the actual data instead.
size_breaks <- function(v) {
  r <- range(v, na.rm = TRUE)
  if (r[1] == r[2]) return(r[1])
  unique(round(seq(r[1], r[2], length.out = 3)))
}

# On Windows a PDF that is currently open in a viewer is locked, and ggsave then
# dies with "cannot open file", killing the whole script. Warn and carry on
# instead - close the viewer and re-run to refresh that one file.
safe_ggsave <- function(file, plot, ...) {
  tryCatch(ggsave(file, plot, ...),
           error = function(e) {
             message("!! could not write ", file, " (", conditionMessage(e), ")")
             message("   -> is it open in a PDF viewer? close it and re-run.")
           })
}

########
# Dotplot
#   y = GO term, ordered by fold enrichment
#   x = fold enrichment (Significant / Expected)
#   point size   = number of candidate genes in the term ("Gene count")
#   point colour = elim Fisher p-value, viridis on a log10 scale
##########
go_dotplot <- function(df, plot_title, file, width = 16, height = NULL, wrap = 45) {
  if (nrow(df) == 0) {
    message("dotplot skipped (no terms): ", file)
    return(invisible(NULL))
  }
  # Height scales with the number of terms, otherwise the y labels overlap.
  if (is.null(height)) height <- max(6, min(60, 0.30 * nrow(df) + 3))

  d <- df %>%
    arrange(FE) %>%
    mutate(Term_wrapped = str_wrap(Term, wrap),
           Term_wrapped = factor(Term_wrapped, levels = unique(Term_wrapped)))

  p <- ggplot(d, aes(x = FE, y = Term_wrapped,
                     size = Significant, colour = elimFisher_num)) +
    geom_point() +
    scale_colour_viridis_c(trans = "log10",
                           labels = function(x) format(x, scientific = TRUE, digits = 1)) +
    scale_size_continuous(breaks = size_breaks(d$Significant), range = c(3, 10)) +
    labs(title = plot_title,
         x = "Fold enrichment (Observed / Expected)", y = NULL,
         colour = "p-value", size = "Gene count") +
    theme_bw() +
    theme(plot.title   = element_text(size = 18, hjust = 0.5),
          axis.text.x  = element_text(size = 16),
          axis.text.y  = element_text(size = 14),
          axis.title.x = element_text(size = 16),
          legend.text  = element_text(size = 13),
          legend.title = element_text(size = 14))
  safe_ggsave(file, p, width = width, height = height, limitsize = FALSE)
  invisible(p)
}

# Bar plot of fold enrichment
go_barplot <- function(df, plot_title, file, width = 16, height = NULL) {
  if (nrow(df) == 0) {
    message("barplot skipped (no terms): ", file)
    return(invisible(NULL))
  }
  if (is.null(height)) height <- max(6, min(60, 0.30 * nrow(df) + 3))
  p <- df %>%
    arrange(FE) %>%    # sorts the data frame but not the factor levels, hence the factor() below
    mutate(Term = factor(Term, levels = Term)) %>%
    ggplot(aes(x = Term, y = FE, fill = elimFisher_num)) +
    geom_bar(position = "dodge", stat = "identity") +
    coord_flip() +
    labs(title = plot_title, y = "Fold enrichment", x = "", fill = "p-value") +
    theme_bw() +
    theme(axis.text.x = element_text(size = 20),
          axis.text.y = element_text(size = 13),
          plot.title = element_text(size = 16),
          axis.title.x = element_text(size = 20),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 20))
  safe_ggsave(file, p, width = width, height = height, limitsize = FALSE)
  invisible(p)
}


########
#3. Run GO enrichment
##########

res <- list()
for (ont in c("BP", "MF", "CC")) {
  r <- run_go(ont)
  a <- add_FE(subset(r$tab, p_num(r$tab$elimFisher) <= 0.05))
  cat(ont, ": feasible genes in universe =", length(genes(r$obj)),
      "| candidates feasible =", numSigGenes(r$obj),
      "| significant terms =", nrow(a), "\n")
  write.table(a, paste0("GO_", ont, "_elim0.05_", tag, ".txt"),
              sep = "\t", quote = FALSE, eol = "\n", row.names = FALSE, col.names = TRUE)
  res[[ont]] <- list(obj = r$obj, sig = a)
}

# For comparison among studies the classic algorithm may be preferable: it does not
# take the GO hierarchy into account and does not discard child terms. Here the
# BP results are ordered by the classic Fisher p-value.
r_cl <- run_go("BP", order_by = "classicFisher")
a_cl <- add_FE(subset(r_cl$tab, p_num(r_cl$tab$elimFisher) <= 0.05))
write.table(a_cl, paste0("GO_BP_classicOrder_", tag, ".txt"),
            sep = "\t", quote = FALSE, eol = "\n", row.names = FALSE, col.names = TRUE)

go_bp <- res$BP$sig
go_mf <- res$MF$sig


################
#4. Plots
################

# merge BP and MF terms
go_all <- rbind(go_bp, go_mf)

# Filter out terms annotated to fewer than 2 or more than MAX_ANNOTATED genes, and terms
# represented by a single candidate gene. There is no single correct way to filter these
# data; adjust the thresholds if needed.
go_all1 <- go_all %>%
  filter(Annotated >= 2) %>%
  filter(Annotated <= MAX_ANNOTATED) %>%
  filter(Significant > 1) %>%
  distinct(Term, .keep_all = TRUE)   # Remove duplicate rows based on Term
cat("BP+MF merged & filtered (2 <= Annotated <=", MAX_ANNOTATED,
    ", Significant > 1):", nrow(go_all1), "terms\n")
write.table(go_all1, paste0("GO_BP_MF_merged_filter_", tag, ".txt"),
            sep = "\t", quote = FALSE, eol = "\n", row.names = FALSE, col.names = TRUE)

go_barplot(go_all1, paste0("GO BP+MF - ", title),
           paste0("GO_BP_MF_merged_filter_barplot_", tag, ".pdf"))
go_dotplot(go_all1, paste0("GO BP+MF significant terms - ", title),
           paste0("GO_BP_MF_merged_filter_dotplot_", tag, ".pdf"))


########
#5. Which of the 61 candidate genes sits under which enriched term
##########
# For every significant BP/MF term: the candidate genes annotated to it, with their
# A. thaliana ortholog and protein name.
genes_in_terms <- function(go.obj, sig_tab, onto) {
  if (nrow(sig_tab) == 0) return(NULL)
  sig <- sigGenes(go.obj)
  ann <- genesInTerm(go.obj, sig_tab$GO.ID)
  do.call(rbind, lapply(seq_len(nrow(sig_tab)), function(i) {
    g <- intersect(ann[[sig_tab$GO.ID[i]]], sig)
    if (length(g) == 0) return(NULL)
    data.frame(ontology = onto, GO.ID = sig_tab$GO.ID[i], Term = sig_tab$Term[i],
               elimFisher = sig_tab$elimFisher[i], gene = g,
               At_ortholog = cand$single_Arabidopsis_thaliana[match(g, cand$gene)],
               protein = cand$UniProt_Protein.names[match(g, cand$gene)],
               stringsAsFactors = FALSE)
  }))
}

term_genes <- rbind(genes_in_terms(res$BP$obj, go_bp, "BP"),
                    genes_in_terms(res$MF$obj, go_mf, "MF"))
if (!is.null(term_genes)) {
  write.table(term_genes, paste0("GO_significant_terms_with_genes_", tag, ".tsv"),
              sep = "\t", quote = FALSE, eol = "\n", row.names = FALSE, col.names = TRUE)
}

sessionInfo()
