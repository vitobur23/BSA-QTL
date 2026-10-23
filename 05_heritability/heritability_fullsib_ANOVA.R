## ---------------------------------------------------------------------------
## Heritability of rosette size in the F2 mapping population, from a full-sib
## (F2-family) variance decomposition (Lynch & Walsh 1998, ch. 18).
##
## F2 individuals are grouped into full-sib families defined by the UNORDERED PAIR of
## their F1 parents. Reciprocal crosses of the same two F1 plants (e.g. family "3030" =
## F1 plant 3030 x plant 3031 and family "3031" = plant 3031 x plant 3030) produce
## genetically equivalent progeny and are therefore merged into one family. The
## `Family` column of the phenotype table is the code of the MOTHER F1 plant; it is
## mapped to the parental pair through the genealogy table.
##
## Quantities reported
##   t       full-sib intraclass correlation, t = sigma2_between / (sigma2_between +
##           sigma2_within). For full sibs t = (1/2 V_A + 1/4 V_D + V_Ec) / V_P, so t is
##           not a heritability: it captures half of the additive and a quarter of the
##           dominance variance and is inflated by the common-family environment V_Ec.
##   2t      the conventional upper bound on narrow-sense heritability (h2 <= 2t),
##           biased upwards by dominance and shared environment.
##   F test  significance of the family effect (one-way ANOVA with the Sokal & Rohlf
##           1995 correction for unequal family sizes).
##   CI      95% bootstrap confidence interval, resampling whole families (5,000 times).
##
## In an F2 design the between-family component captures only the variance due to
## differences among the F1 parental pairs. The Mendelian segregation variance WITHIN
## a family, which is what bulk segregant analysis maps, falls entirely into
## sigma2_within together with environmental variance and measurement error. Hence t is
## a lower bound on the genetic variance available to mapping, and 2t is the relevant
## upper bound.
##
## The script also tests for a reciprocal (maternal / cytoplasmic) effect by nesting
## cross direction within parental pair.
##
## Input files (in ./input; see README):
##   genealogy_table.csv                   F1 parentage and phenotype of the 132 F2
##                                         individuals selected for the bulks
##   F2_phenotypes.csv                     rosette diameter of all phenotyped F2 plants
##                                         (columns Family, Individual, diameter, treatment)
##   per_individual_pairing_with_codes.tsv bulk membership of the sequenced individuals
##
## Analyses
##   B  full, unselected F2 population (n = 762, serpentine soil): PRIMARY result.
##   A  the 132 individuals selected for the bulks. Selecting both tails of the phenotypic
##      distribution inflates the between-family variance, so these estimates are biased
##      upwards and are shown for comparison only.
##   Each analysis is also run with families defined by the directional cross (mother
##   code only), which ignores that reciprocal crosses share both parents.
##
## Run from this folder (setwd() below). Output goes to ./results.
## ---------------------------------------------------------------------------

# setwd("/path/to/05_heritability")

genealogy_path <- file.path("input", "genealogy_table.csv")
pairing_path   <- file.path("input", "per_individual_pairing_with_codes.tsv")
full_pop_path  <- file.path("input", "F2_phenotypes.csv")

out_dir <- "results"
dir.create(out_dir, showWarnings = FALSE)

set.seed(1)
n_boot <- 5000

## ---------------------------------------------------------------------------
## 1. Load data
## ---------------------------------------------------------------------------

geneal_raw <- read.csv(genealogy_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))

pairing <- read.delim(pairing_path, sep = "\t", stringsAsFactors = FALSE)
stopifnot(all(c("bulk", "sequencing_code") %in% names(pairing)))

fullpop_raw <- read.csv(full_pop_path, stringsAsFactors = FALSE, na.strings = c("", "NA"))
stopifnot(all(c("Family", "Individual", "diameter", "treatment") %in% names(fullpop_raw)))

## ---------------------------------------------------------------------------
## 2. Build the analysis tables
## ---------------------------------------------------------------------------

## -- 2a. sample selected for the bulks --
geneal <- geneal_raw[!is.na(geneal_raw[["F2_individual"]]), ]
geneal$f2_id <- as.numeric(geneal[["F2_individual"]])
stopifnot(!any(duplicated(geneal$f2_id)))

geneal$diameter <- as.numeric(geneal[["rosette_diameter"]])
geneal$area     <- as.numeric(geneal[["rosette_area"]])

geneal$mother_code <- as.character(as.numeric(geneal[["F1_mother_code"]]))
geneal$father_code <- as.character(as.numeric(geneal[["F1_father_code"]]))

## directional cross: mother x father
geneal$cross_directional <- paste0("M", geneal$mother_code, "_F", geneal$father_code)
## unordered parental pair = one full-sib family
geneal$pair <- apply(cbind(geneal$mother_code, geneal$father_code), 1,
                     function(x) paste(sort(x), collapse = "x"))

geneal$bulk_corrected <- pairing$bulk[match(geneal$f2_id, pairing$sequencing_code)]

cat("=== (A) Sample selected for the bulks: overview ===\n")
cat(sprintf("F2 individuals with phenotype + parentage: %d\n", nrow(geneal)))
cat(sprintf("  matched to a bulk: %d ; excluded 'middle' individuals: %d\n",
            sum(!is.na(geneal$bulk_corrected)), sum(is.na(geneal$bulk_corrected))))
cat(sprintf("Bulk counts: %s\n",
            paste(names(table(geneal$bulk_corrected)), table(geneal$bulk_corrected),
                  sep = "=", collapse = ", ")))

## ---------------------------------------------------------------------------
## 2b. Crossing design: how many families, and to what extent do they share parents?
## ---------------------------------------------------------------------------

cross_table <- unique(data.frame(mother = geneal$mother_code,
                                 father = geneal$father_code,
                                 mother_name = as.character(geneal[["F1_mother"]]),
                                 father_name = as.character(geneal[["F1_father"]]),
                                 stringsAsFactors = FALSE))
cross_table$pair <- apply(cross_table[, c("mother", "father")], 1,
                          function(x) paste(sort(x), collapse = "x"))
cross_table <- cross_table[order(cross_table$pair, cross_table$mother), ]

n_dir  <- nrow(cross_table)
n_pair <- length(unique(cross_table$pair))
recip_pairs <- names(table(cross_table$pair))[table(cross_table$pair) == 2]

cat("\n=== CROSSING DESIGN ===\n")
cat(sprintf("Directional crosses (mother x father): %d\n", n_dir))
cat(sprintf("Distinct unordered PARENTAL PAIRS  : %d   <-- number of full-sib families\n", n_pair))
cat(sprintf("Pairs crossed in BOTH directions (reciprocal): %d\n", length(recip_pairs)))
cat("Reciprocal pairs:\n"); print(recip_pairs)
cat("\nFull crossing table:\n"); print(cross_table, row.names = FALSE)
write.csv(cross_table, file.path(out_dir, "crossing_design.csv"), row.names = FALSE)

## half-sib links: an F1 plant used as a parent in more than one parental pair
parent_usage <- table(unlist(strsplit(unique(cross_table$pair), "x")))
half_sib_parents <- names(parent_usage)[parent_usage > 1]
cat(sprintf("\nF1 plants in >1 parental pair (=> HALF-SIB links between families): %s\n",
            if (length(half_sib_parents)) paste(half_sib_parents, collapse = ", ") else "none"))
for (p in half_sib_parents) {
  cat(sprintf("  parent %s is in pairs: %s\n", p,
              paste(grep(p, unique(cross_table$pair), value = TRUE), collapse = ", ")))
}

## -- 2c. full, unselected F2 population --
## `Family` there is the mother's code; map it to the unordered parental pair
## via the genealogy. The partner map is built in BOTH directions, so that a
## family recorded only as the reciprocal half (e.g. 3093, whose own row is not
## in the genealogy) still resolves to its pair.
mf <- unique(data.frame(mother = geneal$mother_code, father = geneal$father_code,
                        stringsAsFactors = FALSE))
partner <- c(setNames(mf$father, mf$mother), setNames(mf$mother, mf$father))
partner <- partner[!duplicated(names(partner))]

fullpop <- fullpop_raw[!is.na(fullpop_raw$Individual) & fullpop_raw$Individual != "", ]
fullpop$diameter <- as.numeric(fullpop$diameter)
fullpop$family_mothercode <- as.character(fullpop$Family)   # directional grouping (mother code)
fullpop$pair <- sapply(fullpop$family_mothercode, function(a) {
  if (is.na(partner[a])) NA_character_ else paste(sort(c(a, partner[[a]])), collapse = "x")
})

fullpop_serp    <- fullpop[fullpop$treatment == "serpentine", ]
fullpop_nonserp <- fullpop[fullpop$treatment == "non-serpentine", ]

cat("\n\n=== (B) Full population overview ===\n")
for (nm in c("serpentine", "non-serpentine")) {
  d <- fullpop[fullpop$treatment == nm, ]
  cat(sprintf("%-15s n=%3d ; mother-code groups=%2d ; PARENTAL-PAIR groups=%2d ; unresolved=%d\n",
              nm, nrow(d), length(unique(d$family_mothercode)),
              length(unique(na.omit(d$pair))), sum(is.na(d$pair))))
}
unresolved <- unique(fullpop$family_mothercode[is.na(fullpop$pair)])
if (length(unresolved))
  cat(sprintf("Family codes with no parentage in the genealogy (dropped): %s\n",
              paste(unresolved, collapse = ", ")))

cat("\nParental-pair family sizes (serpentine):\n")
print(sort(table(fullpop_serp$pair)))

## QC: does `Family` really equal the mother's code?
geneal_name <- gsub(" ", "", geneal[["F2_name"]])
serp_name   <- gsub(" ", "", fullpop_serp$Individual)
m <- match(geneal_name, serp_name)
cat(sprintf("\nQC: %d / %d individuals of sample A found in the full serpentine population%s\n",
            sum(!is.na(m)), nrow(geneal),
            if (sum(!is.na(m)) == nrow(geneal)) " -- OK, A is a subset of B." else " -- MISMATCH, investigate."))
disagree <- which(!is.na(m) & fullpop_serp$family_mothercode[m] != geneal$mother_code)
cat(sprintf("QC: `Family` == genealogy mother code for %d / %d individuals.\n",
            sum(!is.na(m)) - length(disagree), sum(!is.na(m))))
if (length(disagree)) {
  cat("    Individuals where the two disagree:\n")
  print(data.frame(individual = geneal_name[disagree],
                   fullpop_Family = fullpop_serp$family_mothercode[m[disagree]],
                   genealogy_mother = geneal$mother_code[disagree],
                   genealogy_father = geneal$father_code[disagree]), row.names = FALSE)
}

## ---------------------------------------------------------------------------
## 3. Full-sib variance-component estimator
##    (Lynch & Walsh 1998 ch.18; Sokal & Rohlf 1995 unequal-n one-way ANOVA)
## ---------------------------------------------------------------------------

fullsib_variance <- function(y, family) {
  ok <- !is.na(y) & !is.na(family)
  y <- y[ok]; family <- factor(family[ok])

  fam_n <- table(family)
  k <- length(fam_n)
  N <- length(y)
  fam_mean <- tapply(y, family, mean)

  ss_between <- sum(fam_n * (fam_mean - mean(y))^2)
  ss_within  <- sum((y - fam_mean[family])^2)
  df_between <- k - 1
  df_within  <- N - k
  ms_between <- ss_between / df_between
  ms_within  <- ss_within / df_within

  n0 <- (N - sum(fam_n^2) / N) / (k - 1)   # Sokal & Rohlf correction, unequal n

  sigma2_within  <- ms_within
  sigma2_between <- max((ms_between - ms_within) / n0, 0)  # variances can't be < 0

  V_P <- sigma2_between + sigma2_within
  t_fs <- if (V_P > 0) sigma2_between / V_P else NA_real_
  f_stat <- ms_between / ms_within

  list(n = N, k_families = k,
       sigma2_between = sigma2_between, sigma2_within = sigma2_within, V_P = V_P,
       t_fullsib_intraclass = t_fs,   # = (1/2 V_A + 1/4 V_D + V_Ec) / V_P
       h2_upper_bound_2t = 2 * t_fs,  # conventional UPPER bound on narrow-sense h2
       F_statistic = f_stat, df_between = df_between, df_within = df_within,
       p_value = pf(f_stat, df_between, df_within, lower.tail = FALSE))
}

## nonparametric bootstrap CI: resample whole families with replacement
bootstrap_ci <- function(y, family, statistic_name, n_boot = 5000, conf = 0.95) {
  ok <- !is.na(y) & !is.na(family)
  y <- y[ok]; family <- factor(family[ok])
  fam_levels <- levels(family)
  idx_by_fam <- lapply(fam_levels, function(f) which(family == f))
  names(idx_by_fam) <- fam_levels

  boot_vals <- numeric(n_boot)
  for (b in seq_len(n_boot)) {
    rs <- sample(fam_levels, length(fam_levels), replace = TRUE)
    idx <- unlist(idx_by_fam[rs])
    # relabel so a family drawn twice counts as two distinct families
    fam_rep <- rep(seq_along(rs), times = lengths(idx_by_fam[rs]))
    res <- tryCatch(fullsib_variance(y[idx], fam_rep), error = function(e) NULL)
    boot_vals[b] <- if (is.null(res)) NA_real_ else res[[statistic_name]]
  }
  alpha <- 1 - conf
  quantile(boot_vals, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
}

run_one <- function(grouping, population_type, sample_name, trait_name, y, family) {
  est <- fullsib_variance(y, family)
  ci_t  <- bootstrap_ci(y, family, "t_fullsib_intraclass", n_boot = n_boot)
  ci_h2 <- bootstrap_ci(y, family, "h2_upper_bound_2t",    n_boot = n_boot)
  data.frame(
    grouping = grouping, population_type = population_type,
    sample = sample_name, trait = trait_name,
    n = est$n, k_families = est$k_families,
    sigma2_between = est$sigma2_between, sigma2_within = est$sigma2_within,
    t_fullsib_intraclass = est$t_fullsib_intraclass,
    t_CI_low = ci_t[1], t_CI_high = ci_t[2],
    h2_upper_bound_2t = est$h2_upper_bound_2t,
    h2_CI_low = ci_h2[1], h2_CI_high = ci_h2[2],
    F_statistic = est$F_statistic, df_between = est$df_between,
    df_within = est$df_within, p_value = est$p_value
  )
}

## ---------------------------------------------------------------------------
## 4. Run all analyses
## ---------------------------------------------------------------------------

results <- list()

## -- (B) full, unselected population: PRIMARY result --
results[["B_serp_pair_diameter"]] <- run_one(
  "parental_pair", "B_full_population_PRIMARY", "serpentine_n762",
  "rosette_diameter", fullpop_serp$diameter, fullpop_serp$pair)

results[["B_serp_directional_diameter"]] <- run_one(
  "mother_code_directional", "B_full_population_PRIMARY", "serpentine_n762",
  "rosette_diameter", fullpop_serp$diameter, fullpop_serp$family_mothercode)

## secondary / G x E context: same families, non-serpentine soil (no shared individuals)
results[["B_nonserp_pair_diameter"]] <- run_one(
  "parental_pair", "B_full_population_secondary_GxE", "non_serpentine",
  "rosette_diameter", fullpop_nonserp$diameter, fullpop_nonserp$pair)

## -- (A) sample selected for the bulks -- upward-biased, comparison only --
trunc_samples <- list(all_available    = geneal,
                      bulked_only_n124 = geneal[!is.na(geneal$bulk_corrected), ])
trunc_traits  <- c(area = "area", rosette_diameter = "diameter")

for (samp_name in names(trunc_samples)) {
  d <- trunc_samples[[samp_name]]
  for (trait_name in names(trunc_traits)) {
    tc <- trunc_traits[[trait_name]]
    results[[paste("A", samp_name, trait_name, "pair", sep = "_")]] <- run_one(
      "parental_pair", "A_truncated_bulk_selected_BIASED",
      samp_name, trait_name, d[[tc]], d$pair)
    results[[paste("A", samp_name, trait_name, "directional", sep = "_")]] <- run_one(
      "mother_code_directional", "A_truncated_bulk_selected_BIASED",
      samp_name, trait_name, d[[tc]], d$cross_directional)
  }
}

results_table <- do.call(rbind, results)
rownames(results_table) <- NULL

cat("\n\n=== Heritability results (full-sib ANOVA, Lynch & Walsh 1998 ch.18) ===\n")
print(results_table[, c("grouping", "population_type", "trait", "n", "k_families",
                        "t_fullsib_intraclass", "t_CI_low", "t_CI_high",
                        "h2_upper_bound_2t", "p_value")], digits = 3)

write.csv(results_table, file.path(out_dir, "heritability_fullsib_results.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4b. Sensitivity of the primary estimate
## ---------------------------------------------------------------------------

sens <- list()
add_sens <- function(label, y, family) {
  e <- fullsib_variance(y, family)
  ci <- bootstrap_ci(y, family, "t_fullsib_intraclass", n_boot = n_boot)
  sens[[label]] <<- data.frame(variant = label, k_families = e$k_families, n = e$n,
                               t_fullsib_intraclass = e$t_fullsib_intraclass,
                               t_CI_low = ci[1], t_CI_high = ci[2],
                               h2_upper_bound_2t = e$h2_upper_bound_2t,
                               p_value = e$p_value)
}
add_sens("primary (all serpentine)", fullpop_serp$diameter, fullpop_serp$pair)
d_no3093 <- fullpop_serp[fullpop_serp$family_mothercode != "3093", ]
add_sens("drop ambiguous family 3093", d_no3093$diameter, d_no3093$pair)
small <- names(table(fullpop_serp$family_mothercode))[table(fullpop_serp$family_mothercode) < 5]
d_nosmall <- fullpop_serp[!fullpop_serp$family_mothercode %in% small, ]
add_sens(sprintf("drop sub-groups n<5 (%s)", paste(small, collapse = ",")),
         d_nosmall$diameter, d_nosmall$pair)
add_sens("sqrt(diameter)", sqrt(fullpop_serp$diameter), fullpop_serp$pair)
add_sens("log(diameter)",  log(fullpop_serp$diameter),  fullpop_serp$pair)

sens_table <- do.call(rbind, sens); rownames(sens_table) <- NULL
cat("\n=== Sensitivity of the primary estimate ===\n")
print(sens_table, digits = 3, row.names = FALSE)
write.csv(sens_table, file.path(out_dir, "heritability_sensitivity.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 4c. Reciprocal (maternal / cytoplasmic) effect
##     = cross direction nested within parental pair
## ---------------------------------------------------------------------------

d <- fullpop_serp[!is.na(fullpop_serp$pair) & !is.na(fullpop_serp$diameter), ]
m_pair <- lm(diameter ~ pair, data = d)
m_dir  <- lm(diameter ~ pair / family_mothercode, data = d)
aov_recip <- anova(m_pair, m_dir)
recip_F  <- aov_recip[["F"]][2]
recip_p  <- aov_recip[["Pr(>F)"]][2]
partial_r2 <- (sum(resid(m_pair)^2) - sum(resid(m_dir)^2)) / sum(resid(m_pair)^2)

cat("\n=== Reciprocal (maternal/cytoplasmic) effect: direction nested within pair ===\n")
print(aov_recip)
cat(sprintf("Partial R2 of cross direction within parental pair: %.4f\n", partial_r2))

recip_means <- aggregate(diameter ~ pair + family_mothercode, data = d,
                         FUN = function(x) c(n = length(x), mean = mean(x)))
recip_means <- do.call(data.frame, recip_means)
names(recip_means)[3:4] <- c("n", "mean_diameter")
recip_means <- recip_means[recip_means$pair %in%
                             names(table(recip_means$pair))[table(recip_means$pair) == 2], ]
recip_means <- recip_means[order(recip_means$pair), ]
cat("\nPer-pair reciprocal means:\n"); print(recip_means, row.names = FALSE, digits = 3)
write.csv(recip_means, file.path(out_dir, "reciprocal_cross_means.csv"), row.names = FALSE)

recip_test <- data.frame(
  test = "cross direction nested within parental pair (serpentine, diameter)",
  n = nrow(d), df_num = aov_recip[["Df"]][2], df_den = aov_recip[["Res.Df"]][2],
  F_statistic = recip_F, p_value = recip_p, partial_R2 = partial_r2)
write.csv(recip_test, file.path(out_dir, "reciprocal_effect_test.csv"), row.names = FALSE)

## ---------------------------------------------------------------------------
## 5. Diagnostic plots
## ---------------------------------------------------------------------------

## 5a. phenotype by parental-pair family, full population
png(file.path(out_dir, "phenotype_by_family_fullpop_diameter.png"),
    width = 1400, height = 700, res = 120)
d2 <- fullpop_serp[!is.na(fullpop_serp$pair), ]
d2$fam <- factor(d2$pair, levels = names(sort(tapply(d2$diameter, d2$pair, mean, na.rm = TRUE))))
par(mar = c(8, 4.5, 3, 1))
boxplot(diameter ~ fam, data = d2, las = 2, cex.axis = 0.75, xlab = "", ylab = "rosette diameter",
        main = sprintf("Rosette diameter by full-sib family (parental pair) -- serpentine, n=%d, k=%d",
                       nrow(d2), length(unique(d2$pair))))
mtext("F2 full-sib family = unordered parental pair (sorted by family mean)", side = 1, line = 6.5)
stripchart(diameter ~ fam, data = d2, vertical = TRUE, method = "jitter",
           pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.25), add = TRUE)
dev.off()

## 5b. phenotype by parental-pair family, sample selected for the bulks (area)
png(file.path(out_dir, "phenotype_by_family_truncated_area.png"),
    width = 1400, height = 700, res = 120)
d1 <- geneal
d1$fam <- factor(d1$pair, levels = names(sort(tapply(d1$area, d1$pair, mean, na.rm = TRUE))))
par(mar = c(8, 4.5, 3, 1))
boxplot(area ~ fam, data = d1, las = 2, cex.axis = 0.75, xlab = "", ylab = "rosette area",
        main = sprintf("Rosette area by full-sib family -- sample selected for the bulks (n=%d, k=%d)",
                       nrow(d1), length(unique(d1$pair))))
mtext("F2 full-sib family = unordered parental pair (sorted by family mean)", side = 1, line = 6.5)
stripchart(area ~ fam, data = d1, vertical = TRUE, method = "jitter",
           pch = 16, cex = 0.5, col = rgb(0, 0, 0, 0.4), add = TRUE)
dev.off()

## 5c. full population distribution with sequenced HIGH/LOW bulk members highlighted
geneal$bulk_label <- ifelse(is.na(geneal$bulk_corrected), "excluded/middle", geneal$bulk_corrected)
bulk_by_name <- setNames(geneal$bulk_label, geneal_name)
fullpop_serp$bulk_membership <- bulk_by_name[serp_name]
fullpop_serp$bulk_membership[is.na(fullpop_serp$bulk_membership)] <- "not sequenced"

png(file.path(out_dir, "fullpop_distribution_with_bulks_highlighted.png"),
    width = 1400, height = 800, res = 120)
cols <- c("not sequenced" = rgb(0.6, 0.6, 0.6, 0.5), "high" = rgb(0.85, 0.1, 0.1, 0.8),
          "low" = rgb(0.1, 0.3, 0.85, 0.8), "excluded/middle" = rgb(0.9, 0.7, 0.1, 0.8))
hist(fullpop_serp$diameter[fullpop_serp$bulk_membership == "not sequenced"],
     breaks = 40, col = cols["not sequenced"], border = NA,
     xlim = range(fullpop_serp$diameter, na.rm = TRUE),
     xlab = "rosette diameter",
     main = "Full F2 population (n=762) with sequenced HIGH/LOW bulks highlighted")
for (grp in c("excluded/middle", "low", "high")) {
  vals <- fullpop_serp$diameter[fullpop_serp$bulk_membership == grp]
  if (length(vals) > 0) rug(vals, col = cols[grp], lwd = 2)
}
legend("topright", legend = names(cols), fill = cols, border = NA, bty = "n")
dev.off()

## 5d. reciprocal-cross means: the maternal/cytoplasmic effect
png(file.path(out_dir, "reciprocal_cross_effect.png"), width = 1200, height = 750, res = 120)
rp <- recip_means
pairs_u <- unique(rp$pair)
ylim <- range(rp$mean_diameter) + c(-0.15, 0.30)
par(mar = c(9, 4.5, 4, 1))
plot(NA, xlim = c(0.5, length(pairs_u) + 0.5), ylim = ylim, xaxt = "n",
     xlab = "", ylab = "family mean rosette diameter",
     main = "Reciprocal cross effect: same two F1 parents, roles swapped")
abline(h = mean(d$diameter), lty = 3, col = "grey50")
for (i in seq_along(pairs_u)) {
  s <- rp[rp$pair == pairs_u[i], ]
  s <- s[order(s$family_mothercode), ]
  segments(i, s$mean_diameter[1], i, s$mean_diameter[2], col = "grey40", lwd = 2)
  points(rep(i, 2), s$mean_diameter, pch = 21, cex = 1.4,
         bg = c(rgb(0.85, 0.35, 0.1, 0.9), rgb(0.1, 0.4, 0.75, 0.9)))
  text(rep(i, 2), s$mean_diameter, labels = s$family_mothercode, pos = 4, cex = 0.55)
}
axis(1, at = seq_along(pairs_u), labels = pairs_u, las = 2, cex.axis = 0.8)
mtext("parental pair (each point = one cross direction, labelled by mother's code)",
      side = 1, line = 7.5)
legend("topleft", bty = "n", cex = 0.85,
       legend = sprintf("direction within pair: F = %.2f, p = %.2g, partial R2 = %.3f",
                        recip_F, recip_p, partial_r2))
dev.off()

cat(sprintf("\nResults and plots written to: %s\n", out_dir))

## ---------------------------------------------------------------------------
## 6. Summary of the key numbers
## ---------------------------------------------------------------------------

pri <- results[["B_serp_pair_diameter"]]
dirn <- results[["B_serp_directional_diameter"]]
cat("\n=== Key numbers ===\n")
cat(sprintf("PRIMARY (parental pair): t = %.3f (95%% CI %.3f-%.3f), 2t = %.3f, k = %d, n = %d, p = %.2g\n",
            pri$t_fullsib_intraclass, pri$t_CI_low, pri$t_CI_high,
            pri$h2_upper_bound_2t, pri$k_families, pri$n, pri$p_value))
cat(sprintf("Directional grouping (mother code): t = %.3f, k = %d  -- %.0f%% higher, because the\n  reciprocal-cross effect is counted as between-family variance\n",
            dirn$t_fullsib_intraclass, dirn$k_families,
            100 * (dirn$t_fullsib_intraclass / pri$t_fullsib_intraclass - 1)))
cat("\n- t is the full-sib INTRACLASS CORRELATION, not H2. For full sibs\n")
cat("  t = (1/2 V_A + 1/4 V_D + V_Ec)/V_P, so 2t is the upper bound on narrow-sense heritability.\n")
cat("- In an F2 design the within-family Mendelian segregation variance (what BSA maps)\n")
cat("  sits in sigma2_within, so t is a LOWER bound on the genetic variance available to mapping.\n")
cat("- The A rows (sample selected for the bulks) are upward-biased by the selection of the two\n")
cat("  phenotypic tails and are shown for comparison only.\n")
cat(sprintf("- Design: %d directional crosses among %d parental pairs, %d of them reciprocal;\n  half-sib links via F1 plant(s) %s.\n",
            n_dir, n_pair, length(recip_pairs), paste(half_sib_parents, collapse = ", ")))
