# =============================================================================
# 04 — Score A: discrimination, DeLong, Firth model, bootstrap internal validation,
#      calibration, thresholds, sensitivity populations, site analyses
# =============================================================================
dat <- readRDS(file.path(OUT_DIR, "wegio_analysis_dataset.rds"))


# ---- 1. AUC of Score A and Beales count; DeLong -----------------------------------------------------------
auc_A      <- auc_boot_ci(dat$scoreA, dat$bbs)
auc_Bcount <- auc_boot_ci(dat$beales_count, dat$bbs)
rocA <- roc(dat$bbs, dat$scoreA, direction = "<", quiet = TRUE)
rocB <- roc(dat$bbs, dat$beales_count, direction = "<", quiet = TRUE)
dl   <- roc.test(rocA, rocB, method = "delong", paired = TRUE)
dl_diff <- unname(auc(rocA) - auc(rocB)); dl_se <- dl_diff / unname(dl$statistic)
set.seed(SEED_BOOT); i1 <- which(dat$bbs == 1); i0 <- which(dat$bbs == 0)
diff_bt <- replicate(N_BOOT_AUC, { s <- c(rs(i1), rs(i0))
  auc_fast(dat$scoreA[s], dat$bbs[s]) - auc_fast(dat$beales_count[s], dat$bbs[s]) })
delong_res <- tibble(
  auc_scoreA = fmt_auc(auc_A), auc_scoreA_delongCI = paste(sprintf("%.3f", ci.auc(rocA, method = "delong")[c(1, 3)]), collapse = "–"),
  auc_beales_count = fmt_auc(auc_Bcount), diff = dl_diff,
  diff_delong_CI = sprintf("%.3f–%.3f", dl_diff - 1.96 * dl_se, dl_diff + 1.96 * dl_se),
  diff_boot_CI = sprintf("%.3f–%.3f", quantile(diff_bt, .025), quantile(diff_bt, .975)),
  delong_z = unname(dl$statistic), delong_p = dl$p.value)

# ---- 2. Firth model (seven items + age + male) ---------------------------------------------------------------------
dat$age10 <- dat$age_y / 10
dat$age10[is.na(dat$age10)] <- median(dat$age10, na.rm = TRUE)     # one missing age: median imputation
FORM  <- bbs ~ poly + retinal + hypog_genital + renal + synbrachy + nd + dental + age10 + male
ITEMS <- c("poly", "retinal", "hypog_genital", "renal", "synbrachy", "nd", "dental")

fit_firth <- function(d, pl = TRUE, flic = FALSE)
  logistf(FORM, data = d, pl = pl, flic = flic, control = logistf.control(maxit = 200), plcontrol = logistf.control(maxit = 200))
derive_points <- function(beta) ifelse(beta > 0.25, pmax(1, round(beta)), 0)   # manuscript rule
lp_of <- function(fit, d) drop(model.matrix(FORM, d) %*% coef(fit))

fit <- fit_firth(dat)
fit_flic <- fit_firth(dat, pl = FALSE, flic = TRUE)
b  <- coef(fit); se <- sqrt(diag(vcov(fit)))
table6_model <- tibble(term = names(b), beta = b, OR = exp(b),
  OR_wald_lo = exp(b - 1.96 * se), OR_wald_hi = exp(b + 1.96 * se),
  OR_pl_lo = exp(fit$ci.lower), OR_pl_hi = exp(fit$ci.upper), p_pl = fit$prob) %>%
  mutate(scoreA_points = c(NA, SCORE_A_POINTS[c("poly", "retinal", "hypog_genital", "renal", "synbrachy", "nd")], 0, NA, NA),
         points_reestimated = ifelse(term %in% ITEMS, derive_points(beta), NA))
pts0 <- setNames(table6_model$points_reestimated[match(ITEMS, table6_model$term)], ITEMS)
score_from_points <- function(d, pts) drop(as.matrix(d[, ITEMS]) %*% pts[ITEMS])
dat$score_rederived <- score_from_points(dat, pts0)

# ---- 3. Bootstrap internal validation (Harrell optimism, 500 resamples) ------------------------------------
brier <- function(p, y) mean((p - y)^2)
cal_slope <- function(lp, y) unname(coef(glm(y ~ lp, family = binomial))[2])
cal_int   <- function(lp, y) unname(coef(glm(y ~ offset(lp), family = binomial))[1])
lp_app <- lp_of(fit, dat); p_app <- plogis(lp_app)
app <- c(auc_model = auc_fast(lp_app, dat$bbs), auc_score = auc_fast(dat$score_rederived, dat$bbs),
         slope = 1, brier = brier(p_app, dat$bbs), citl = cal_int(lp_app, dat$bbs))
set.seed(SEED_BOOT)
boot <- map_dfr(seq_len(N_BOOT), function(i) {
  s <- sample(nrow(dat), replace = TRUE); db <- dat[s, ]
  fb <- tryCatch(fit_firth(db, pl = FALSE), error = function(e) NULL); if (is.null(fb)) return(NULL)
  pts <- derive_points(coef(fb)[ITEMS])
  lp_b <- lp_of(fb, db); lp_o <- lp_of(fb, dat)
  sc_b <- score_from_points(db, pts); sc_o <- score_from_points(dat, pts)
  tibble(b = i,
         opt_auc_model = auc_fast(lp_b, db$bbs) - auc_fast(lp_o, dat$bbs),
         opt_auc_score = auc_fast(sc_b, db$bbs) - auc_fast(sc_o, dat$bbs),
         slope_test = cal_slope(lp_o, dat$bbs),
         opt_brier = brier(plogis(lp_b), db$bbs) - brier(plogis(lp_o), dat$bbs),
         !!!setNames(as.list(pts), paste0("pts_", ITEMS)))
})
validation <- tibble(
  n_boot_ok = nrow(boot),
  auc_model_apparent = app["auc_model"], auc_model_corrected = app["auc_model"] - mean(boot$opt_auc_model),
  auc_score_rederived_apparent = app["auc_score"], auc_score_rederived_corrected = app["auc_score"] - mean(boot$opt_auc_score),
  auc_scoreA_fixed_apparent = auc_A[1],
  cal_slope_corrected = mean(boot$slope_test),
  brier_apparent = app["brier"], brier_corrected = app["brier"] - mean(boot$opt_brier),
  citl_firth_apparent = app["citl"],
  citl_flic_apparent = cal_int(lp_of(fit_flic, dat), dat$bbs),
  mean_pred_firth = mean(p_app), mean_pred_flic = mean(plogis(lp_of(fit_flic, dat))), observed = mean(dat$bbs))
tableS14 <- map_dfr(ITEMS, function(it) { v <- boot[[paste0("pts_", it)]]
  tibble(item = it, scoreA = c(SCORE_A_POINTS, dental = 0)[it], pct_ge1 = 100 * mean(v >= 1),
         median = median(v), p2.5 = quantile(v, .025), p97.5 = quantile(v, .975),
         pct_equal_scoreA = 100 * mean(v == c(SCORE_A_POINTS, dental = 0)[it])) })

# ---- 4. Calibration (deciles + smoothed) -----------------------------------------------------------------------------
calib <- dat %>% mutate(p = p_app, p_flic = plogis(lp_of(fit_flic, dat)),
                        dec = ntile(p, 10)) %>% group_by(dec) %>%
  summarise(n = n(), pred = mean(p), pred_flic = mean(p_flic), obs = mean(bbs), events = sum(bbs),
            lo = binom.wilson(sum(bbs), n())$lower, hi = binom.wilson(sum(bbs), n())$upper)
saveRDS(list(calib = calib, p = p_app, p_flic = plogis(lp_of(fit_flic, dat)), y = dat$bbs), file.path(OUT_DIR, "calibration.rds"))

# ---- 5. Thresholds (table 6 lower) and Score A bands (table 7) --------------------------------------------------------
thr_tab <- function(d, score = "scoreA", ks = 1:13) map_dfr(ks, function(k) {
  t <- as.integer(d[[score]] >= k); tp <- sum(t & d$bbs == 1); fp <- sum(t & d$bbs == 0)
  fn <- sum(!t & d$bbs == 1); tn <- sum(!t & d$bbs == 0)
  tibble(threshold = k, n = tp + fp, pct = 100 * (tp + fp) / nrow(d), bbs = tp,
         sens = fmt_ci(tp, tp + fn), spec = fmt_ci(tn, tn + fp),
         ppv = if (tp + fp > 0) fmt_ci(tp, tp + fp) else NA, npv = round(100 * tn / (tn + fn), 1)) })
table6_thr <- thr_tab(dat)
table7 <- dat %>% mutate(band = cut(scoreA, c(-Inf, 2, 5, 8, Inf), labels = c("0–2", "3–5", "6–8", "≥9"))) %>%
  group_by(band) %>% summarise(n = n(), bbs = sum(bbs), nonBBS_ob = sum(grp4 == 2), nonob = sum(grp4 == 3), none = sum(grp4 == 4),
    p_bbs = fmt_ci(bbs, n), p_obgd = fmt_ci(bbs + nonBBS_ob, n), p_gd = fmt_ci(bbs + nonBBS_ob + nonob, n))

# ---- 6. Table 8 / S7: reference standards and sensitivity populations -------------------------------------------------
perf_pop <- function(d, label) {
  if (sum(d$bbs) == 0) return(NULL)
  a <- auc_boot_ci(d$scoreA, d$bbs)
  f <- function(v) { tp <- sum(d[[v]] & d$bbs == 1); sp <- 100 * mean(d[[v]][d$bbs == 0] == 0)
    sprintf("%d (%.1f); %.1f", tp, 100 * tp / sum(d$bbs), sp) }
  tibble(population = label, N = nrow(d), BBS = sum(d$bbs), auc = fmt_auc(a),
         scoreA3 = f("scoreA3"), scoreA5 = f("scoreA5"), beales = f("beales"), ern = f("ern"))
}
catA <- dat %>% mutate(bbs = as.integer(bbs_cat %in% "A"))
table8 <- bind_rows(
  perf_pop(dat, "All study-assigned BBS"),
  perf_pop(filter(catA, is.na(bbs_cat) | bbs_cat == "A"), "Category A; other BBS excluded"),
  perf_pop(dat %>% filter(is.na(bbs_cat) | bbs_cat %in% c("A", "B1")), "A+B1; below-rule excluded"),
  perf_pop(catA, "Category A vs all others"),
  perf_pop(filter(dat, index == 1), "Index participants only"),
  perf_pop(filter(dat, high_site == 0), "Excluding sites 10 and 39"),
  perf_pop(filter(dat, high_site == 1), "Sites 10 and 39 only"),
  perf_pop(filter(dat, age18 %in% "<18"), "<18 years"),
  perf_pop(filter(dat, age18 %in% ">=18"), ">=18 years"),
  perf_pop(filter(dat, age16 %in% "<16"), "<16 years"),
  perf_pop(filter(dat, age16 %in% ">=16"), ">=16 years"),
  perf_pop(filter(dat, clin_dx == 0), "No clinical BBS diagnosis before enrolment"),
  perf_pop(filter(dat, pd == 0), "No protocol deviation"),
  perf_pop(filter(dat, age_y < 6), "<6 years"))
s7_extra <- tibble(
  analysis = c("Retinal item removed, all", "Retinal item removed, <16", "Retinal item removed, <18", "Retinal item removed, >=18",
               "Not-assessed counted as present", "Complete case (all items assessed)"),
  N = c(1004, sum(dat$age16 %in% "<16"), sum(dat$age18 %in% "<18"), sum(dat$age18 %in% ">=18"), 1004, sum(dat$scoreA_complete)),
  BBS = c(38, sum(dat$bbs[dat$age16 %in% "<16"]), sum(dat$bbs[dat$age18 %in% "<18"]), sum(dat$bbs[dat$age18 %in% ">=18"]), 38,
          sum(dat$bbs[dat$scoreA_complete])),
  auc = c(auc_fast(dat$scoreA_noret, dat$bbs),
          with(filter(dat, age16 %in% "<16"), auc_fast(scoreA_noret, bbs)),
          with(filter(dat, age18 %in% "<18"), auc_fast(scoreA_noret, bbs)),
          with(filter(dat, age18 %in% ">=18"), auc_fast(scoreA_noret, bbs)),
          auc_fast(dat$scoreA_napres, dat$bbs),
          with(filter(dat, scoreA_complete), auc_fast(scoreA, bbs))),
  sens3 = c(mean(dat$scoreA_noret[dat$bbs == 1] >= 3), NA, NA, NA, mean(dat$scoreA_napres[dat$bbs == 1] >= 3), NA) * 100,
  spec3 = c(mean(dat$scoreA_noret[dat$bbs == 0] < 3), NA, NA, NA, mean(dat$scoreA_napres[dat$bbs == 0] < 3), NA) * 100)
retinal_spec <- dat %>% filter(!is.na(age18), bbs == 0) %>% group_by(age18) %>% summarise(spec_retinal = 100 * mean(retinal == 0))

# ---- 7. Site analyses -------------------------------------------------------------------------------------------------
loso <- map_dfr(sort(unique(dat$site)), function(s) { d <- filter(dat, site != s)
  tibble(site_left_out = s, n = nrow(d), bbs = sum(d$bbs), auc = auc_fast(d$scoreA, d$bbs)) })
loso_range <- range(loso$auc)

# ---- 8. Predictive values at lower prevalence (table S15) ---------------------------------------------------------------
se5 <- mean(dat$scoreA5[dat$bbs == 1]); sp5 <- mean(dat$scoreA5[dat$bbs == 0] == 0)
tableS15 <- tibble(prev = c(.001, .005, .01, .02, mean(dat$bbs), .05)) %>%
  mutate(ppv = 100 * se5 * prev / (se5 * prev + (1 - sp5) * (1 - prev)),
         npv = 100 * sp5 * (1 - prev) / (sp5 * (1 - prev) + (1 - se5) * prev))

# ---- 9. Thresholds by age (table S8) --------------------------------------------------------------------------------------
tableS8 <- bind_rows(thr_tab(filter(dat, age16 %in% "<16"), ks = c(2:6, 8)) %>% mutate(stratum = "<16"),
                     thr_tab(filter(dat, age16 %in% ">=16"), ks = c(2:6, 8)) %>% mutate(stratum = ">=16"),
                     thr_tab(filter(dat, age18 %in% "<18"), ks = c(2:6, 8)) %>% mutate(stratum = "<18"),
                     thr_tab(filter(dat, age18 %in% ">=18"), ks = c(2:6, 8)) %>% mutate(stratum = ">=18"))

saveRDS(list(rocA = rocA, rocB = rocB, dat = dat, fit = fit, boot = boot), file.path(OUT_DIR, "scoreA_objects.rds"))
write_xlsx(list(delong = delong_res, table6_model = table6_model, table6_thresholds = table6_thr, table7 = table7,
                validation = validation, tableS14 = tableS14, calibration_deciles = calib, table8 = table8,
                tableS7_extra = s7_extra, retinal_specificity = retinal_spec, loso = loso, tableS15 = tableS15,
                tableS8 = tableS8, bootstrap_raw = boot),
           file.path(OUT_DIR, "tables", "04_scoreA_model.xlsx"))
message("04 Score A done; LOSO AUC range ", paste(round(loso_range, 3), collapse = "–"))
