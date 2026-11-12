# =============================================================================
# 10 — Cell-by-cell comparison of every rebuilt table (09) with the tables in the manuscript .docx
#      (accepted view of tracked changes). Set MANUSCRIPT_DOCX to the file to check.
#      Output: output/tables/10_docx_vs_R_mismatches.xlsx
# =============================================================================
source(file.path(DATA_DIR, "R_workflow_V25", "R", "read_docx_tables.R"))
MANUSCRIPT_DOCX <- Sys.getenv("MANUSCRIPT_DOCX", unset = file.path(DATA_DIR, "WEGIO_manuscript_20261112_V25_2_TDTC.docx"))
DT <- read_docx_tables(MANUSCRIPT_DOCX)
MT <- readRDS(file.path(OUT_DIR, "manuscript_tables.rds"))

# manuscript table -> position of the table in the .docx (0-based, document order)
MAP <- c(table1 = 0, table2 = 1, table4 = 4, table5 = 5, table6_upper = 6, table6_lower = 7, table7 = 8, table8 = 9,
         tableS2 = 11, tableS4 = 12, tableS5 = 13, tableS6_upper = 14, tableS6_lower = 15, tableS7 = 16, tableS8 = 17,
         tableS11 = 20, tableS13 = 22, tableS14 = 23, tableS15 = 24, tableS16 = 25)
norm <- function(s) gsub("\\s+", " ", trimws(gsub("−", "-", s)))

cmp <- map_dfr(names(MAP), function(nm) {
  doc_rows <- DT[[MAP[[nm]] + 1]][-1]                            # drop header row
  doc_rows <- doc_rows[vapply(doc_rows, function(r) any(nzchar(r)), logical(1))]
  R <- MT[[nm]]
  if (nm == "tableS4") {                                          # match participants by ID
    ids <- vapply(doc_rows, `[`, "", 1); ord <- match(R[, 1], ids)
    if (anyNA(ord)) warning("S4: IDs not found in docx")
    doc_rows <- doc_rows[ord]
  }
  n <- max(length(doc_rows), nrow(R))
  map_dfr(seq_len(n), function(i) {
    dr <- if (i <= length(doc_rows)) doc_rows[[i]] else character(0)
    rr <- if (i <= nrow(R)) R[i, ] else character(0)
    k <- max(length(dr), length(rr))
    tibble(table = nm, row = i, col = seq_len(k),
           docx = norm(c(dr, rep(NA, k - length(dr)))), R = norm(c(rr, rep(NA, k - length(rr)))))
  })
}) %>% mutate(status = case_when(is.na(R) ~ "input (not rebuilt)", docx == R ~ "match", TRUE ~ "MISMATCH"))

summary_tab <- cmp %>% count(table, status) %>% tidyr::pivot_wider(names_from = status, values_from = n, values_fill = 0)
mism <- filter(cmp, status == "MISMATCH")
write_xlsx(list(summary = summary_tab, mismatches = mism, all_cells = cmp), file.path(OUT_DIR, "tables", "10_docx_vs_R_mismatches.xlsx"))
message("10 compare: ", nrow(cmp), " cells checked; ", nrow(mism), " mismatches (", basename(MANUSCRIPT_DOCX), ")")
