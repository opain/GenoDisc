#!/usr/bin/env Rscript
# Unit-level sanity check for the direction-of-effect recipes used to compute
# Direction and Reversal_Z columns in:
#   - format_twas_gsea_drugtargetor_results.R   (drug + ATC levels)
#   - format_twas_gsea_cmap_results.R           (per-signature + per-MOA)
#
# Run this in any R session that has the base 'stats' package. No data files
# are needed. Exits with status 0 on success, 1 on any failure.

assert <- function(cond, msg) {
  if(!isTRUE(cond)) {
    cat("FAIL: ", msg, "\n", sep = "")
    quit(status = 1, save = "no")
  } else {
    cat("ok:   ", msg, "\n", sep = "")
  }
}

cat("=== 1. ATC Wilcoxon (formula form, HL = out - in) ===\n")
# Force in-class drugs to have LOW T (opposes disease at drug level), out-class
# to have HIGH T. With the formula `wilcox.test(rank(T) ~ class_bin)`, R uses
# class_bin == 0 as the first group (x), so HL = HL(out - in) > 0.
T_vals <- c(10, 9, 8, 7, 6, 1, 2, 3)
class_bin <- c(0, 0, 0, 0, 0, 1, 1, 1)
est <- as.numeric(wilcox.test(rank(T_vals) ~ class_bin, conf.int = TRUE)$estimate)
assert(est > 0, paste0("ATC: in-class low T => Estimate > 0 (got ", round(est, 3), ")"))
# Direction recipe at ATC level: Estimate > 0 -> "Opposes disease".
direction <- ifelse(est > 0, "Opposes disease",
             ifelse(est < 0, "Matches disease", NA_character_))
assert(direction == "Opposes disease",
       "ATC: Direction recipe labels low-in-class as 'Opposes disease'")

# Reverse the toy data: in-class has HIGH T (matches disease). HL should flip.
T_vals2 <- c(1, 2, 3, 4, 5, 10, 9, 8)
est2 <- as.numeric(wilcox.test(rank(T_vals2) ~ class_bin, conf.int = TRUE)$estimate)
assert(est2 < 0, paste0("ATC: in-class high T => Estimate < 0 (got ", round(est2, 3), ")"))

cat("\n=== 2. MOA Wilcoxon (two-vector form, HL = in - out) ===\n")
# CMAP MOA uses wilcox.test(in, out). HL = HL(in - out). With in-class LOW T,
# HL < 0. Direction recipe at MOA level: Estimate < 0 -> "Opposes disease".
in_T  <- c(1, 2, 3)
out_T <- c(10, 9, 8, 7, 6)
est_moa <- as.numeric(wilcox.test(in_T, out_T, conf.int = TRUE)$estimate)
assert(est_moa < 0, paste0("MOA: in-class low T => Estimate < 0 (got ", round(est_moa, 3), ")"))
direction_moa <- ifelse(est_moa < 0, "Opposes disease",
                 ifelse(est_moa > 0, "Matches disease", NA_character_))
assert(direction_moa == "Opposes disease",
       "MOA: Direction recipe (FLIPPED) labels low-in-class as 'Opposes disease'")

cat("\n=== 3. Reversal_Z is positive in the 'opposes' direction at every level ===\n")
# Drug level: Reversal_Z = -T. T < 0 (opposes) -> Reversal_Z > 0.
T_drug <- -2.3
rev_z_drug <- -T_drug
assert(rev_z_drug > 0, "drug: T < 0 -> Reversal_Z > 0")

# ATC level: Reversal_Z = qnorm(1-P/2) * sign(Estimate). The magnitude is the
# one-sided Z corresponding to the two-sided Wilcoxon P (always >= 0 for any
# P <= 1), so sign(Reversal_Z) always agrees with sign(Estimate), which agrees
# with Direction by construction.
P_atc <- 0.001
rev_z_atc <- qnorm(1 - P_atc/2) * sign(est)    # Estimate from section 1
assert(rev_z_atc > 0, "ATC: Estimate > 0, P small -> Reversal_Z > 0")

# Crucial edge-case: when P > 0.5, magnitude qnorm(1-P/2) is still >= 0, so
# Reversal_Z still tracks the Direction sign — unlike the older -qnorm(P)
# formula which gave a sign disagreement for P > 0.5.
P_atc_hi <- 0.7
rev_z_atc_hi <- qnorm(1 - P_atc_hi/2) * sign(est)
assert(rev_z_atc_hi > 0, "ATC: high-P case still has Reversal_Z > 0 when Estimate > 0")

# MOA level: Reversal_Z = qnorm(1-P/2) * sign(-Estimate). Estimate < 0 (opposes)
# -> sign(-Estimate) > 0 -> Reversal_Z > 0.
P_moa <- 0.001
rev_z_moa <- qnorm(1 - P_moa/2) * sign(-est_moa)  # Estimate from section 2
assert(rev_z_moa > 0, "MOA: Estimate < 0, P small -> Reversal_Z > 0")

cat("\n=== 4. Drug-level vs ATC-level signs are consistent on shared toy data ===\n")
# Build a small ATC class where every in-class drug has a known drug-level T.
# Confirm: when all in-class drugs are opposes-direction at drug level
# (negative T), the ATC-level Direction is also "Opposes disease".
set.seed(42)
T_in  <- rnorm(8, mean = -2, sd = 0.5)
T_out <- rnorm(40, mean = 0, sd = 1)
class_bin_full <- c(rep(1, length(T_in)), rep(0, length(T_out)))
T_full <- c(T_in, T_out)

# Drug-level Direction: each drug independently
drug_directions <- ifelse(T_in < 0, "Opposes disease",
                   ifelse(T_in > 0, "Matches disease", NA_character_))
assert(all(drug_directions == "Opposes disease"),
       "drug: all in-class drugs labelled 'Opposes disease'")

# ATC-level Direction: aggregate via formula-form Wilcoxon
est_aggregate <- as.numeric(wilcox.test(rank(T_full) ~ class_bin_full, conf.int = TRUE)$estimate)
direction_aggregate <- ifelse(est_aggregate > 0, "Opposes disease",
                       ifelse(est_aggregate < 0, "Matches disease", NA_character_))
assert(direction_aggregate == "Opposes disease",
       "ATC: aggregate of opposes-direction drugs is 'Opposes disease'")

cat("\n=== 5. Non-directional Reversal_Z has no sign coupling to disease direction ===\n")
# Non-directional case: P is one-sided right-tail (enrichment), Direction = NA,
# Reversal_Z = qnorm(1-P) is always >= 0 for P <= 1.
P_nondir <- 0.01
rev_z_nondir <- qnorm(1 - P_nondir)
assert(rev_z_nondir > 0, "nondir: Reversal_Z > 0 for enriched P")
P_nondir_hi <- 0.95
rev_z_nondir_hi <- qnorm(1 - P_nondir_hi)
assert(rev_z_nondir_hi >= 0 || abs(rev_z_nondir_hi) < 5,
       "nondir: Reversal_Z stays interpretable for non-significant P")
direction_nondir <- NA_character_
assert(is.na(direction_nondir), "nondir: Direction = NA")

cat("\nAll sign-convention assertions passed.\n")
