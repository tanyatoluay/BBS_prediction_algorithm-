# =============================================================================
# 04b — Validation of the complete Score A derivation procedure
#   Runs after 04_scoreA_model.R in the same session (uses dat, FORM, ITEMS, fit_firth,
#   derive_points, score_from_points).
#   Procedure being validated (as in the original derivation, Wegio-Biostats.Rmd):
#     prespecified 7 clinical items -> Firth model (+ age, sex) -> integer points
#     (round(beta), >=1 if beta > 0.25) -> score -> threshold
#   A. Bootstrap optimism (500) for AUC AND for sensitivity/specificity at the thresholds
#      (fixed >=3 / >=5, and a threshold re-chosen in every resample by Youden's index)
#   B. Repeated stratified 10-fold cross-validation (20 repeats) of the whole procedure
#   C. Internal-external validation by site: derive at sites 10/39, test at all other sites
# =============================================================================

THR <- c(THR_LO, THR_HI)
sens_spec <- function(score, y, k) c(sens = mean(score[y == 1] >= k), spec = mean(score[y == 0] < k))
youden_k  <- function(score, y) { ks <- sort(unique(score)); ks <- ks[ks > 0]
  j <- sapply(ks, function(k) sum(sens_spec(score, y, k)) - 1); ks[which.max(j)] }
derive_on <- function(d) { fb <- fit_firth(d, pl = FALSE); derive_points(coef(fb)[ITEMS]) }

# apparent values: fixed Score A and the score re-derived on the full cohort
app_tab <- map_dfr(c(scoreA = "scoreA", legacy = "scoreA_legacy"), function(v) {
  s <- dat[[v]]; k_y <- youden_k(s, dat$bbs)
  tibble(score = v, auc = auc_fast(s, dat$bbs),
         sens3 = sens_spec(s, dat$bbs, THR_LO)["sens"], spec3 = sens_spec(s, dat$bbs, THR_LO)["spec"],
         senshi = sens_spec(s, dat$bbs, THR_HI)["sens"], spechi = sens_spec(s, dat$bbs, THR_HI)["spec"],
         youden_threshold = k_y, sensY = sens_spec(s, dat$bbs, k_y)["sens"], specY = sens_spec(s, dat$bbs, k_y)["spec"])
})

# ---- A. bootstrap optimism of the full procedure -------------------------------------------------
set.seed(SEED_BOOT)
bootB <- map_dfr(seq_len(N_BOOT), function(i) {
  db <- dat[sample(nrow(dat), replace = TRUE), ]
  pts <- tryCatch(derive_on(db), error = function(e) NULL); if (is.null(pts)) return(NULL)
  sb <- score_from_points(db, pts); so <- score_from_points(dat, pts)
  ky <- youden_k(sb, db$bbs)
  ssb <- function(k) sens_spec(sb, db$bbs, k); sso <- function(k) sens_spec(so, dat$bbs, k)
  tibble(opt_auc = auc_fast(sb, db$bbs) - auc_fast(so, dat$bbs),
         opt_sens3 = ssb(THR_LO)["sens"] - sso(THR_LO)["sens"], opt_spec3 = ssb(THR_LO)["spec"] - sso(THR_LO)["spec"],
         opt_senshi = ssb(THR_HI)["sens"] - sso(THR_HI)["sens"], opt_spechi = ssb(THR_HI)["spec"] - sso(THR_HI)["spec"],
         youden_k = ky, opt_sensY = ssb(ky)["sens"] - sso(ky)["sens"], opt_specY = ssb(ky)["spec"] - sso(ky)["spec"])
})
opt <- colMeans(select(bootB, starts_with("opt_")))
corrected <- app_tab %>% mutate(
  auc_corr = auc - opt["opt_auc"],
  sens3_corr = sens3 - opt["opt_sens3"], spec3_corr = spec3 - opt["opt_spec3"],
  senshi_corr = senshi - opt["opt_senshi"], spechi_corr = spechi - opt["opt_spechi"],
  sensY_corr = sensY - opt["opt_sensY"], specY_corr = specY - opt["opt_specY"])
youden_dist <- count(bootB, youden_k) %>% mutate(pct = 100 * n / sum(n))

# ---- B. repeated stratified 10-fold cross-validation -------------------------------------------------
set.seed(SEED_BOOT)
cv <- map_dfr(1:20, function(r) {
  fold <- integer(nrow(dat))
  for (cl in 0:1) { idx <- which(dat$bbs == cl); fold[idx] <- sample(rep_len(1:10, length(idx))) }
  pred <- numeric(nrow(dat))
  for (f in 1:10) { pts <- derive_on(dat[fold != f, ]); pred[fold == f] <- score_from_points(dat[fold == f, ], pts) }
  tibble(rep = r, auc = auc_fast(pred, dat$bbs),
         sens3 = sens_spec(pred, dat$bbs, THR_LO)["sens"], spec3 = sens_spec(pred, dat$bbs, THR_LO)["spec"],
         senshi = sens_spec(pred, dat$bbs, THR_HI)["sens"], spechi = sens_spec(pred, dat$bbs, THR_HI)["spec"])
})
cv_summary <- cv %>% summarise(across(-rep, list(mean = mean, min = min, max = max)))

# ---- C. internal-external validation by site ------------------------------------------------------------
dev <- filter(dat, high_site == 1); val <- filter(dat, high_site == 0)
pts_site <- derive_on(dev); sv <- score_from_points(val, pts_site)
site_ext <- tibble(
  derivation = "sites 10 + 39", n_dev = nrow(dev), bbs_dev = sum(dev$bbs),
  validation = "35 other sites", n_val = nrow(val), bbs_val = sum(val$bbs),
  points = paste(names(pts_site), pts_site, sep = "=", collapse = ", "),
  auc_val = auc_fast(sv, val$bbs), auc_val_ci = fmt_auc(auc_boot_ci(sv, val$bbs)),
  senshi = sprintf("%d/%d", sum(sv[val$bbs == 1] >= THR_HI), sum(val$bbs)),
  spechi = 100 * mean(sv[val$bbs == 0] < THR_HI),
  auc_fixed_scoreA_val = auc_fast(val$scoreA, val$bbs))

write_xlsx(list(apparent_and_corrected = corrected, youden_threshold_distribution = youden_dist,
                cv_repeats = cv, cv_summary = cv_summary, site_internal_external = site_ext,
                bootstrap_procedure_raw = bootB),
           file.path(OUT_DIR, "tables", "04b_full_validation.xlsx"))
message("04b full-procedure validation done")
