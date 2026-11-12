# =============================================================================
# WEGIO V25 reproducible analysis workflow — configuration
# Inputs : study genetic database V21.0 (VSK, 4 Nov 2026) + eCRF export V4.0 (PBE, 1 Nov 2026)
# Author : T. D. Toluay Coşkun (workflow), generated 2026-11-09
# =============================================================================

invisible(Sys.setlocale("LC_ALL", "en_US.UTF-8"))
suppressPackageStartupMessages({
  library(readxl); library(dplyr); library(tidyr); library(stringr); library(purrr)
  library(logistf); library(pROC); library(binom); library(exact2x2)
  library(ggplot2); library(patchwork); library(writexl)
})

# ---- paths -----------------------------------------------------------------
# Run from the project folder ("WEGIO new grouping") or set WEGIO_DIR.
DATA_DIR <- Sys.getenv("WEGIO_DIR", unset = getwd())
GEN_FILE  <- file.path(DATA_DIR, "WEGIO_all_genetic_results_20261104_V21.0_VSK.xlsx")
ECRF_FILE <- file.path(DATA_DIR, "eCRF_Wegio_Export_Sheet_20261101_V4.0_PBE.xlsx")
OUT_DIR   <- file.path(DATA_DIR, "R_workflow_V25", "output")
dir.create(file.path(OUT_DIR, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "tables"),  recursive = TRUE, showWarnings = FALSE)

# ---- analysis constants ------------------------------------------------------
SEED_BOOT   <- 20261104   # manuscript table 6 footnote (table S14 states 20261105 -> harmonise)
N_BOOT      <- 500        # Firth internal validation
N_BOOT_AUC  <- 2000       # bootstrap CI for AUCs (table 8)
N_PERM      <- 5000       # Monte Carlo exact test
HIGH_BBS_SITES <- c("10", "39")

# Score A — re-derived on V21.0 / eCRF V4.0 (Firth, round(beta), >=1 point if beta > 0.25; 04c)
SCORE_A_POINTS  <- c(poly = 5, retinal = 2, hypog_genital = 1, renal = 1, synbrachy = 1, nd = 1)   # max 11
# legacy Score A (earlier analysis, used in manuscript V25) — kept for the V25 reproduction check
SCORE_A_LEGACY  <- c(poly = 5, retinal = 2, hypog_genital = 2, renal = 2, synbrachy = 1, nd = 1)   # max 13
THR_LO <- 3   # highest threshold retaining maximal sensitivity (04c)
THR_HI <- 4   # maximum optimism-corrected Youden index (04c)

# ---- BBS evidence categories (manuscript table S11, as specified by VSK 27 Oct 2026) ----
# Categories are not a column in V21.0. They are read from a LOCAL file kept next to the data
# (never committed): bbs_evidence_categories.csv with columns id, category (A, B1, B2, C).
# Template: config/bbs_evidence_categories_TEMPLATE.csv. Cross-checked against V21.0 VUS tiers in 01.
CAT_FILE <- file.path(DATA_DIR, "bbs_evidence_categories.csv")
if (!file.exists(CAT_FILE)) stop("Missing ", CAT_FILE, " (see config/bbs_evidence_categories_TEMPLATE.csv)")
cat_tab <- read.csv(CAT_FILE, colClasses = "character")
BBS_CAT <- split(trimws(cat_tab$id), factor(cat_tab$category, levels = c("A", "B1", "B2", "C")))

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- helpers -------------------------------------------------------------------
wilson <- function(x, n) {
  if (n == 0) return(c(est = NA, lo = NA, hi = NA))
  r <- binom::binom.wilson(x, n)
  c(est = 100 * x / n, lo = 100 * r$lower, hi = 100 * r$upper)
}
fmt_ci <- function(x, n, d = 1) {
  w <- wilson(x, n); sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), w[1], w[2], w[3])
}
npct <- function(x, n, d = 1) sprintf(paste0("%d (%.", d, "f%%)"), x, 100 * x / n)
med_iqr <- function(v, d = 1) {
  q <- quantile(v, c(.5, .25, .75), na.rm = TRUE, type = 7)
  sprintf(paste0("%.", d, "f (%.", d, "f–%.", d, "f)"), q[1], q[2], q[3])
}

# ---- fast AUC (Mann–Whitney) and stratified bootstrap CI ---------------------------------------
auc_fast <- function(score, y) {
  y <- as.integer(y); n1 <- sum(y == 1); n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  r <- rank(score); (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}
rs <- function(x) x[sample.int(length(x), length(x), replace = TRUE)]   # safe resample
auc_boot_ci <- function(score, y, B = N_BOOT_AUC, seed = SEED_BOOT) {
  set.seed(seed); i1 <- which(y == 1); i0 <- which(y == 0)
  bt <- replicate(B, { s <- c(rs(i1), rs(i0)); auc_fast(score[s], y[s]) })
  c(auc = auc_fast(score, y), lo = unname(quantile(bt, .025)), hi = unname(quantile(bt, .975)))
}
fmt_auc <- function(a) sprintf("%.3f (%.3f–%.3f)", a[1], a[2], a[3])
# number with mid-dot decimal separator (manuscript style)
md <- function(x, d = 1) ifelse(is.na(x), "—", chartr(".", "·", formatC(x, format = "f", digits = d)))
