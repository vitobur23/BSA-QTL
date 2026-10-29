############################
# Tissue-specific expression of the candidate genes
#
# Expression of the candidate genes (via their A. thaliana orthologs) in root tip, root
# and leaf tissues, from the Arabidopsis tissue atlas (Expression Atlas E-MTAB-7978).
# For each tissue, genes of the whole atlas are classified as lowly (<= 25th percentile
# of mean TPM), medium or highly (>= 75th percentile) expressed; the proportions of the
# candidate genes in these classes are then compared with the genome-wide proportions
# expected for a random set of genes.
#
# Run from the folder that contains this script, after putting the input files into ./input
# (see README). Output files are written to ./results.
############################

# setwd("/path/to/07_gene_expression")

library(readr)
library(data.table)
library(dplyr)
library(tidyr)

############################
# File paths
############################

# Expression Atlas download (TPM table), see the link below
atlas_file <- file.path("input", "E-MTAB-7978-query-results.tpms.tsv")

# Candidate genes: tab-separated table with a column of A. thaliana gene identifiers
candidate_file <- file.path("input", "Candidate_genes_AT_codes.txt")

output_directory <- "results"

# Column containing A. thaliana gene identifiers
candidate_column <- "A.thaliana"

dir.create(
  output_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


############################
# Import expression atlas
############################

# Input file available from:
# https://www.ebi.ac.uk/gxa/experiments/E-MTAB-7978/Downloads

tissue_expression_atlas <- data.table(
  read_tsv(
    atlas_file,
    comment = "#",
    skip = 4,
    show_col_types = FALSE
  )
)

# Rename the first column
setnames(
  tissue_expression_atlas,
  old = names(tissue_expression_atlas)[1],
  new = "Gene_ID"
)

# Standardize gene identifiers
tissue_expression_atlas[
  ,
  Gene_ID := toupper(trimws(as.character(Gene_ID)))
]


############################
# Identify tissue types
############################

expression_columns <- setdiff(
  names(tissue_expression_atlas),
  "Gene_ID"
)

tissue_types <- data.table(
  variable = expression_columns,
  tissue = trimws(
    vapply(
      strsplit(expression_columns, ", "),
      function(x) {
        if (length(x) >= 2) x[2] else x[1]
      },
      character(1)
    )
  )
)

# Inspect the tissue labels
print(sort(unique(tissue_types$tissue)))


############################
# Transform to long format
############################

tissue_expression_long <- melt(
  tissue_expression_atlas,
  id.vars = "Gene_ID",
  measure.vars = expression_columns,
  variable.name = "variable",
  value.name = "TPM"
)

tissue_expression_long <- merge(
  tissue_expression_long,
  tissue_types,
  by = "variable",
  all.x = TRUE
)

# Convert TPM values to numeric and replace missing values with zero
tissue_expression_long[
  ,
  TPM := as.numeric(TPM)
]

tissue_expression_long[
  is.na(TPM),
  TPM := 0
]


############################
# Calculate mean expression
############################

# Mean TPM is first calculated per gene and original tissue type.
tissue_expression_sum <- tissue_expression_long[
  ,
  .(
    tpm_mean = mean(TPM, na.rm = TRUE),
    tpm_max = max(TPM, na.rm = TRUE)
  ),
  by = .(Gene_ID, tissue)
]

write.table(
  tissue_expression_sum,
  file.path(
    output_directory,
    "Tissue_expression_mean_E-MTAB-7978.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################
# Select analyzed tissues
############################

leaf_tissues <- c(
  "rosette leaf",
  "cauline leaf",
  "adult vascular leaf"
)

# Root
tissue_exp_root <- tissue_expression_sum[
  tissue == "root",
  .(Gene_ID, tpm_mean, tpm_max)
]

# Root tip
tissue_exp_root_tip <- tissue_expression_sum[
  tissue == "root tip",
  .(Gene_ID, tpm_mean, tpm_max)
]

# Leaves
#
# The three leaf tissue types are combined into one leaf-expression
# value per gene to avoid counting each gene multiple times.
tissue_exp_leaf <- tissue_expression_sum[
  tissue %in% leaf_tissues,
  .(
    tpm_mean = mean(tpm_mean, na.rm = TRUE),
    tpm_max = max(tpm_max, na.rm = TRUE),
    n_leaf_tissues = uniqueN(tissue)
  ),
  by = Gene_ID
]

# Verify the selected tissues
cat("\nNumbers of genes in the selected tissue categories:\n")
cat("Root:", nrow(tissue_exp_root), "\n")
cat("Root tip:", nrow(tissue_exp_root_tip), "\n")
cat("Leaf:", nrow(tissue_exp_leaf), "\n")


############################
# Classify expression levels
############################

classify_expression <- function(expression_data, tissue_name) {
  
  low_threshold <- quantile(
    expression_data$tpm_mean,
    probs = 0.25,
    na.rm = TRUE,
    names = FALSE
  )
  
  high_threshold <- quantile(
    expression_data$tpm_mean,
    probs = 0.75,
    na.rm = TRUE,
    names = FALSE
  )
  
  expression_data %>%
    mutate(
      Tissue = tissue_name,
      low_threshold = low_threshold,
      high_threshold = high_threshold,
      Expression_level = case_when(
        tpm_mean <= low_threshold ~ "Low",
        tpm_mean >= high_threshold ~ "High",
        TRUE ~ "Medium"
      ),
      Expression_level = factor(
        Expression_level,
        levels = c("Low", "Medium", "High")
      )
    )
}


############################
# Apply expression classes
############################

tissue_exp_root <- classify_expression(
  tissue_exp_root,
  "Root"
)

tissue_exp_root_tip <- classify_expression(
  tissue_exp_root_tip,
  "Root tip"
)

tissue_exp_leaf <- classify_expression(
  tissue_exp_leaf,
  "Leaf"
)

# Combine genome-wide classifications
all_tissue_expression <- bind_rows(
  tissue_exp_root_tip,
  tissue_exp_root,
  tissue_exp_leaf
)


############################
# Import new candidate genes
############################

candidate_raw <- read.table(
  candidate_file,
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (!candidate_column %in% names(candidate_raw)) {
  stop(
    paste0(
      "Column '", candidate_column,
      "' was not found in the candidate file. ",
      "Available columns: ",
      paste(names(candidate_raw), collapse = ", ")
    )
  )
}

# Extract canonical A. thaliana identifiers, e.g. AT1G01010.
# This also removes transcript suffixes such as AT1G01010.1.
candidate_genes <- candidate_raw %>%
  transmute(
    Gene_ID = toupper(
      sub(
        ".*(AT[1-5CM]G[0-9]{5}).*",
        "\\1",
        as.character(.data[[candidate_column]])
      )
    )
  ) %>%
  filter(
    grepl(
      "^AT[1-5CM]G[0-9]{5}$",
      Gene_ID
    )
  ) %>%
  distinct(Gene_ID)

cat(
  "\nNumber of unique candidate genes:",
  nrow(candidate_genes),
  "\n"
)


############################
# Merge candidates with expression
############################

candidate_expression <- all_tissue_expression %>%
  inner_join(
    candidate_genes,
    by = "Gene_ID"
  )

cat("\nCandidate genes matched in each tissue:\n")

print(
  candidate_expression %>%
    group_by(Tissue) %>%
    summarise(
      n_candidates = n_distinct(Gene_ID),
      .groups = "drop"
    )
)


############################
# Genome-wide null expectations
############################

# These genome-wide proportions represent the proportions expected
# if a random set of genes were selected from the expression atlas.
null_proportions <- all_tissue_expression %>%
  count(
    Tissue,
    Expression_level,
    name = "n_all_genes"
  ) %>%
  complete(
    Tissue,
    Expression_level = factor(
      c("Low", "Medium", "High"),
      levels = c("Low", "Medium", "High")
    ),
    fill = list(n_all_genes = 0)
  ) %>%
  group_by(Tissue) %>%
  mutate(
    n_all_genes_total = sum(n_all_genes),
    expected_proportion = n_all_genes / n_all_genes_total,
    expected_percentage = 100 * expected_proportion
  ) %>%
  ungroup()


############################
# Observed candidate proportions
############################

candidate_proportions <- candidate_expression %>%
  count(
    Tissue,
    Expression_level,
    name = "n_candidates"
  ) %>%
  complete(
    Tissue,
    Expression_level = factor(
      c("Low", "Medium", "High"),
      levels = c("Low", "Medium", "High")
    ),
    fill = list(n_candidates = 0)
  ) %>%
  group_by(Tissue) %>%
  mutate(
    n_candidates_total = sum(n_candidates),
    candidate_proportion =
      n_candidates / n_candidates_total,
    candidate_percentage =
      100 * candidate_proportion
  ) %>%
  ungroup()


############################
# Combine observed and expected proportions
############################

expression_summary <- candidate_proportions %>%
  left_join(
    null_proportions,
    by = c("Tissue", "Expression_level")
  ) %>%
  mutate(
    expected_candidate_count =
      n_candidates_total * expected_proportion,
    
    proportion_difference =
      candidate_proportion - expected_proportion,
    
    fold_relative_to_random =
      candidate_proportion / expected_proportion
  ) %>%
  transmute(
    Tissue,
    Expression_level = as.character(Expression_level),
    n_candidates,
    n_candidates_total,
    candidate_percentage =
      round(candidate_percentage, 1),
    n_all_genes,
    n_all_genes_total,
    expected_random_percentage =
      round(expected_percentage, 1),
    expected_candidate_count =
      round(expected_candidate_count, 1),
    fold_relative_to_random =
      round(fold_relative_to_random, 2)
  ) %>%
  arrange(
    factor(
      Tissue,
      levels = c("Root tip", "Root", "Leaf")
    ),
    factor(
      Expression_level,
      levels = c("High", "Medium", "Low")
    )
  )

cat("\nCandidate proportions and random-set expectations:\n")

print(expression_summary)


############################
# Write principal output files
############################

# Gene-level candidate classifications
write.table(
  candidate_expression,
  file.path(
    output_directory,
    "Candidate_tissue_expression_classes.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Summary containing candidate and null proportions
write.table(
  expression_summary,
  file.path(
    output_directory,
    "Candidate_expression_and_null_expectations.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################
# Write original tissue outputs
############################

candidate_root_tip <- candidate_expression %>%
  filter(Tissue == "Root tip")

candidate_root <- candidate_expression %>%
  filter(Tissue == "Root")

candidate_leaf <- candidate_expression %>%
  filter(Tissue == "Leaf")

write.table(
  candidate_root_tip,
  file.path(
    output_directory,
    "Tissue_expression_E-MTAB-7978_root_tip.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  candidate_root,
  file.path(
    output_directory,
    "Tissue_expression_E-MTAB-7978_root.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

write.table(
  candidate_leaf,
  file.path(
    output_directory,
    "Tissue_expression_E-MTAB-7978_leaves.txt"
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)


############################
# Print manuscript values
############################

cat("\nHighly expressed genes:\n")

print(
  expression_summary %>%
    filter(Expression_level == "High") %>%
    select(
      Tissue,
      n_candidates,
      n_candidates_total,
      candidate_percentage,
      expected_random_percentage
    )
)

cat("\nLowly expressed genes:\n")

print(
  expression_summary %>%
    filter(Expression_level == "Low") %>%
    select(
      Tissue,
      n_candidates,
      n_candidates_total,
      candidate_percentage,
      expected_random_percentage
    )
)