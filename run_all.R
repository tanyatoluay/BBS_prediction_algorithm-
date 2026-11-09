# =============================================================================
# WEGIO V25 — run the complete reproducible workflow
#   Rscript R_workflow_V25/run_all.R          (from the "WEGIO new grouping" folder)
# Outputs: R_workflow_V25/output/{tables,figures}, sessionInfo.txt
# =============================================================================
if (!file.exists("WEGIO_all_genetic_results_20261104_V21.0_VSK.xlsx"))
  stop("Run from the folder that contains the V21.0 and eCRF V4.0 workbooks (or set WEGIO_DIR).")
need <- c("readxl", "dplyr", "tidyr", "stringr", "purrr", "logistf", "pROC", "binom", "exact2x2",
          "ggplot2", "patchwork", "writexl")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")

R_DIR <- file.path("R_workflow_V25", "R")
for (f in c("00_config.R", "01_build_dataset.R", "02_descriptives.R", "03_criteria.R", "04_scoreA_model.R",
            "05_age_route_nonBBS.R", "06_figures.R", "07_check_manuscript.R")) {
  message("== ", f); source(file.path(R_DIR, f), echo = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
