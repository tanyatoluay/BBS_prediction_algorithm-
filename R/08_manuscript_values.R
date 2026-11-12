# =============================================================================
# 08 — All Score A–dependent values for the manuscript (re-derived Score A, thresholds THR_LO/THR_HI)
#      Runs after 01–05 in the same session. Writes output/manuscript_values.json (formatted, mid-dot decimals)
# =============================================================================
md <- function(x, d = 1) chartr(".", "·", formatC(x, format = "f", digits = d))
pc <- function(k, n, d = 1) md(100 * k / n, d)
ci3 <- function(a) sprintf("%s (%s–%s)", md(a[1], 3), md(a[2], 3), md(a[3], 3))
wil <- function(k, n) { w <- wilson(k, n); sprintf("%s (%s–%s)", md(w[1]), md(w[2]), md(w[3])) }
fp <- function(p) if (p < 0.0001) "<0·0001" else if (p < 0.001) md(p, 4) else if (p < 0.01) md(p, 3) else md(p, 2)
d <- dat; y <- d$bbs; s <- d$scoreA; nb <- sum(y == 0)
k_lo <- THR_LO; k_hi <- THR_HI
V <- list()
acc <- function(sc, yy, k) { tp <- sum(sc[yy == 1] >= k); fn <- sum(yy == 1) - tp; fp_ <- sum(sc[yy == 0] >= k); tn <- sum(yy == 0) - fp_
  list(tp = tp, n1 = tp + fn, sens = wil(tp, tp + fn), spec = wil(tn, tn + fp_), sel = tp + fp_,
       ppv = if (tp + fp_ > 0) wil(tp, tp + fp_) else "—", npv = md(100 * tn / (tn + fn)), sens_p = pc(tp, tp + fn), spec_p = pc(tn, tn + fp_),
       ppv_p = if (tp + fp_ > 0) pc(tp, tp + fp_) else "—") }

# --- discrimination, DeLong ---------------------------------------------------------------------
V$auc <- ci3(auc_A); V$auc_pt <- md(auc_A[1], 3)
V$auc_beales <- md(auc_Bcount[1], 3)
V$diff <- md(delong_res$diff, 3); V$diff_ci <- sub("–", "–", chartr(".", "·", delong_res$diff_boot_CI)); V$delong_p <- fp(delong_res$delong_p)
V$auc_corr <- md(validation$auc_score_rederived_corrected, 3); V$auc_model_app <- md(validation$auc_model_apparent, 3)
V$auc_model_corr <- md(validation$auc_model_corrected, 3); V$slope <- md(validation$cal_slope_corrected, 2)
V$brier <- md(validation$brier_apparent, 3); V$brier_corr <- md(validation$brier_corrected, 3)
# --- thresholds ------------------------------------------------------------------------------------
a_lo <- acc(s, y, k_lo); a_hi <- acc(s, y, k_hi)
V$lo <- a_lo; V$hi <- a_hi
V$lo_sel_pct <- pc(a_lo$sel, nrow(d)); V$hi_sel_pct <- pc(a_hi$sel, nrow(d))
V$below_lo <- sprintf("%d of %d", sum(y[s < k_lo]), sum(s < k_lo)); V$below_lo_pct <- pc(sum(y[s < k_lo]), sum(s < k_lo))
se <- mean(s[y == 1] >= k_hi); sp <- mean(s[y == 0] < k_hi)
ppv_at <- function(p) 100 * se * p / (se * p + (1 - sp) * (1 - p)); npv_at <- function(p) 100 * sp * (1 - p) / (sp * (1 - p) + (1 - se) * p)
V$ppv_prev <- setNames(lapply(c(.001, .005, .01, .02, mean(y), .05), function(p) c(prev = md(100 * p), ppv = md(ppv_at(p)), npv = md(npv_at(p)))),
                       c("0.1", "0.5", "1", "2", "obs", "5"))
bb <- filter(d, bbs == 1)
V$mcn_only_score <- sum(bb$scoreAhi == 1 & bb$beales == 0); V$mcn_only_beales <- sum(bb$scoreAhi == 0 & bb$beales == 1)
V$mcn_p <- fp(binom.test(V$mcn_only_score, V$mcn_only_score + V$mcn_only_beales)$p.value)
V$spec_beales <- pc(sum(d$beales[y == 0] == 0), nb)
V$polydactyly_alone <- sum(d$poly == 1 & d$scoreA >= k_hi) == sum(d$poly)    # 5 points alone reach the threshold
V$n_poly <- sum(d$poly); V$n_poly_bbs <- sum(d$poly[y == 1])
# --- category A / populations ------------------------------------------------------------------
gA <- d %>% filter(is.na(bbs_cat) | bbs_cat == "A") %>% mutate(bbs = as.integer(bbs_cat %in% "A"))
V$catA_auc <- ci3(auc_boot_ci(gA$scoreA, gA$bbs)); V$catA_hi <- sum(gA$scoreA[gA$bbs == 1] >= k_hi)
V$index_auc <- md(with(filter(d, index == 1), auc_fast(scoreA, bbs)), 3)
V$loso <- sprintf("%s to %s", md(min(loso$auc), 3), md(max(loso$auc), 3))
oth <- filter(d, high_site == 0); V$oth_auc <- ci3(auc_boot_ci(oth$scoreA, oth$bbs)); V$oth_hi <- sum(oth$scoreA[oth$bbs == 1] >= k_hi)
V$lt18 <- acc(d$scoreA[d$age18 %in% "<18"], y[d$age18 %in% "<18"], k_hi); V$ge18 <- acc(d$scoreA[d$age18 %in% ">=18"], y[d$age18 %in% ">=18"], k_hi)
V$noret_adult_auc <- md(with(filter(d, age18 %in% ">=18"), auc_fast(scoreA_noret, bbs)), 3)
V$napres_auc <- md(auc_fast(d$scoreA_napres, y), 3)
V$cc_n <- sum(d$scoreA_complete); V$cc_bbs <- sum(y[d$scoreA_complete]); V$cc_auc <- md(with(filter(d, scoreA_complete), auc_fast(scoreA, bbs)), 3)
# --- full-procedure validation (04b) ---------------------------------------------------------------
cr <- corrected %>% filter(score == "scoreA")
V$proc_auc <- md(cr$auc_corr, 3); V$proc_sens_hi <- md(100 * cr$senshi_corr); V$proc_spec_hi <- md(100 * cr$spechi_corr)
V$proc_sens_lo <- md(100 * cr$sens3_corr); V$proc_spec_lo <- md(100 * cr$spec3_corr)
V$proc_sensY <- md(100 * cr$sensY_corr); V$proc_specY <- md(100 * cr$specY_corr)
V$youden_pct_hi <- md(youden_dist$pct[youden_dist$youden_k == k_hi], 0)
V$youden_pct_3to5 <- md(sum(youden_dist$pct[youden_dist$youden_k %in% 3:5]), 0)
V$cv_auc <- md(cv_summary$auc_mean, 3); V$cv_range <- sprintf("%s–%s", md(cv_summary$auc_min, 3), md(cv_summary$auc_max, 3))
V$cv_sens_hi <- md(100 * cv_summary$senshi_mean); V$cv_spec_hi <- md(100 * cv_summary$spechi_mean)
V$site_auc <- md(site_ext$auc_val, 3); V$site_ci <- sub("^[^(]*\\(", "", sub("\\)$", "", chartr(".", "·", site_ext$auc_val_ci)))
V$site_hi <- site_ext$senshi; V$site_spec <- md(site_ext$spechi)
# site-cluster bootstrap (04d)
cc_ <- function(st, pct = FALSE, d = 3) { r <- cluster_ci[cluster_ci$statistic == st, ]
  if (pct) sprintf("%s–%s", md(100 * r$cluster_lo, 1), md(100 * r$cluster_hi, 1)) else sprintf("%s–%s", md(r$cluster_lo, d), md(r$cluster_hi, d)) }
V$cl_auc <- cc_("auc"); V$cl_diff <- cc_("diff"); V$cl_sens_hi <- cc_("sens_hi", TRUE); V$cl_spec_hi <- cc_("spec_hi", TRUE)
V$cl_sens_lo <- cc_("sens_lo", TRUE); V$cl_spec_lo <- cc_("spec_lo", TRUE); V$n_sites <- length(unique(d$site))
# threshold-choice table (04c)
tc <- threshold_choice %>% filter(score == "rederived_V21")
V$thr_youden <- setNames(lapply(seq_len(nrow(tc)), function(i) c(k = tc$threshold[i], sens = md(100 * tc$sens_corr[i]), spec = md(100 * tc$spec_corr[i]), J = md(tc$youden_corr[i], 3))), tc$threshold)

# --- tables -----------------------------------------------------------------------------------------------
# Table 4 rows (>=lo, >=hi) x strata
T4 <- table4 %>% filter(criterion %in% c(paste0("Score A ≥", k_lo), paste0("Score A ≥", k_hi))) %>%
  transmute(criterion, stratum, identified, sens = chartr(".", "·", sens), spec = chartr(".", "·", spec), selected, ppv = md(ppv), npv = md(npv))
# Table 6 lower: thresholds 1..max
T6 <- map_dfr(1:sum(SCORE_A_POINTS), function(k) { a <- acc(s, y, k)
  tibble(k = paste0("≥", k), n = a$sel, pct = pc(a$sel, nrow(d)), bbs = a$tp, sens = a$sens, spec = a$spec, ppv = a$ppv, npv = a$npv) })
# Table 7 bands
T7 <- table7 %>% transmute(band = as.character(band), n, bbs, p_bbs = chartr(".", "·", p_bbs), nonBBS_ob, p_obgd = chartr(".", "·", p_obgd),
                           nonob, p_gd = chartr(".", "·", p_gd), none)
# Table 8
pop8 <- list(all = d, catA_excl = filter(gA, TRUE), ab1 = d %>% filter(is.na(bbs_cat) | bbs_cat %in% c("A", "B1")),
             catA_all = d %>% mutate(bbs = as.integer(bbs_cat %in% "A")), index = filter(d, index == 1), oth = oth,
             lt18 = filter(d, age18 %in% "<18"), ge18 = filter(d, age18 %in% ">=18"))
T8 <- imap_dfr(pop8, function(p, nm) { lo <- acc(p$scoreA, p$bbs, k_lo); hi <- acc(p$scoreA, p$bbs, k_hi)
  tibble(pop = nm, auc = ci3(auc_boot_ci(p$scoreA, p$bbs)), lo = sprintf("%d (%s)", lo$tp, lo$sens_p),
         hi = sprintf("%d (%s); %s", hi$tp, hi$sens_p, hi$spec_p)) })
# S4: score per BBS participant
S4 <- bb %>% transmute(id, scoreA)
# S7 rows (>= THR_LO)
s7row <- function(p, sc = "scoreA") { a <- acc(p[[sc]], p$bbs, k_lo); c(N = nrow(p), BBS = sum(p$bbs), auc = md(auc_fast(p[[sc]], p$bbs), 3),
  sens = a$sens_p, spec = a$spec_p, ppv = a$ppv_p, npv = a$npv) }
mol_ex <- d %>% filter(!(bbs_cat %in% c("B2", "C"))); mol_nc <- d %>% mutate(bbs = as.integer(bbs_cat %in% c("A", "B1")))
S7 <- list(primary = s7row(d), mol_excl = s7row(mol_ex), mol_noncase = s7row(mol_nc), no_cdx = s7row(filter(d, clin_dx == 0)),
           no_pd = s7row(filter(d, pd == 0)), lt16 = s7row(filter(d, age16 %in% "<16")),
           lt16_noret = s7row(filter(d, age16 %in% "<16"), "scoreA_noret"), ge16 = s7row(filter(d, age16 %in% ">=16")),
           noret = s7row(d, "scoreA_noret"), lt6 = s7row(filter(d, age_y < 6)), napres = s7row(d, "scoreA_napres"),
           cc = s7row(filter(d, scoreA_complete)))
# S8 thresholds by age16 (2,3,4,5,6,8)
S8 <- map_dfr(c("<16", ">=16"), function(a) { p <- filter(d, age16 %in% a); map_dfr(c(2:6, 8), function(k) { x <- acc(p$scoreA, p$bbs, k)
  tibble(stratum = a, k = paste0("≥", k), n = x$sel, bbs = x$tp, sens = x$sens_p, spec = x$spec_p, ppv = x$ppv_p, npv = x$npv) }) })
V$s8_auc_lt16 <- md(with(filter(d, age16 %in% "<16"), auc_fast(scoreA, bbs)), 3); V$s8_auc_ge16 <- md(with(filter(d, age16 %in% ">=16"), auc_fast(scoreA, bbs)), 3)
# S11 Score A >= hi by category ; by molecular status
S11 <- bb %>% group_by(bbs_cat) %>% summarise(hi = sum(scoreAhi), lo = sum(scoreA3))
V$mol_rule_lo <- sum(bb$scoreA3[bb$bbs_cat %in% c("A", "B1")]); V$below_rule_lo <- sum(bb$scoreA3[bb$bbs_cat %in% c("B2", "C")])
# S14
S14 <- tableS14 %>% transmute(item, pts = scoreA, ge1 = md(pct_ge1), med = sprintf("%d (%d–%d)", median, as.integer(p2.5), as.integer(p97.5)), eq = md(pct_equal_scoreA))
# Table 6 upper
T6u <- table6_model %>% filter(term %in% ITEMS) %>% transmute(term, new = SCORE_A_POINTS[term], legacy = SCORE_A_LEGACY[term])

jsonlite::write_json(list(V = V, T4 = T4, T6 = T6, T6u = T6u, T7 = T7, T8 = T8, S4 = S4, S7 = S7, S8 = S8, S11 = S11, S14 = S14),
                     file.path(OUT_DIR, "manuscript_values.json"), auto_unbox = TRUE, pretty = TRUE, na = "string")
message("08 manuscript values written")
