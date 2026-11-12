# =============================================================================
# WEGIO V25 — run the complete reproducible workflow
#   Rscript R_workflow_V25/run_all.R          (from the "WEGIO new grouping" folder)
# Outputs: R_workflow_V25/output/{tables,figures}, sessionInfo.txt
# =============================================================================
if (!file.exists("WEGIO_all_genetic_results_20261104_V21.0_VSK.xlsx"))
  stop("Run from the folder that contains the V21.0 and eCRF V4.0 workbooks (or set WEGIO_DIR).")
need <- c("readxl", "dplyr", "tidyr", "stringr", "purrr", "logistf", "pROC", "binom", "exact2x2",
          "ggplot2", "patchwork", "writexl", "jsonlite")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")

R_DIR <- file.path("R_workflow_V25", "R")
for (f in c("00_config.R", "01_build_dataset.R", "02_descriptives.R", "03_criteria.R", "04_scoreA_model.R", "04b_full_validation.R", "04c_rederive_score.R", "04d_cluster_bootstrap.R",
            "05_age_route_nonBBS.R", "06_figures.R", "07_results_values.R", "08_results_tables.R")) {
  message("== ", f); source(file.path(R_DIR, f), echo = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
# README figures (aggregate results only) are refreshed from the current run
fig_docs <- c("Figure2_clinical_profile.png", "Figure3_scoreA_criteria.png", "FigureS1_calibration_manuscript.png", "FigureS2_adults_by_route.png")
dir.create(file.path("R_workflow_V25", "docs", "figures"), recursive = TRUE, showWarnings = FALSE)
invisible(file.copy(file.path(OUT_DIR, "figures", fig_docs), file.path("R_workflow_V25", "docs", "figures"), overwrite = TRUE))
