# =============================================================================
# 04c — Score selection on the current data (V21.0 + eCRF V4.0)
#   Runs after 04 / 04b in the same session.
#   Decision 1 (points): integer points re-derived from the Firth model on the current data
#                        (rule: round(beta), >=1 point if beta > 0.25) vs legacy Score A points.
#   Decision 2 (threshold): for every threshold 2..9, bootstrap optimism-corrected sensitivity and
#                        specificity of the COMPLETE procedure (points re-derived in each resample,
#                        threshold fixed); choose the threshold with the highest corrected Youden index.
#   Decision 3 (stability): how often each item keeps its point value / >=1 point across resamples.
# =============================================================================

KS <- 2:9
SCORES <- c(legacy_scoreA = "scoreA_legacy", rederived_V21 = "scoreA")
pts_v21 <- pts0                                   # from 04: points re-derived on current data
dat$score_v21 <- dat$score_rederived

ss_k <- function(s, y, k) c(sens = mean(s[y == 1] >= k), spec = mean(s[y == 0] < k))

# apparent performance per threshold
apparent_k <- map_dfr(names(SCORES), function(nm) map_dfr(KS, function(k) {
  s <- dat[[SCORES[[nm]]]]; v <- ss_k(s, dat$bbs, k)
  tibble(score = nm, threshold = k, selected = sum(s >= k), bbs_identified = sum(s[dat$bbs == 1] >= k),
         sens_app = v["sens"], spec_app = v["spec"], ppv_app = mean(dat$bbs[s >= k]))
}))

# optimism of the complete procedure for each threshold
set.seed(SEED_BOOT)
opt_k <- map_dfr(seq_len(N_BOOT), function(i) {
  db <- dat[sample(nrow(dat), replace = TRUE), ]
  pts <- tryCatch(derive_on(db), error = function(e) NULL); if (is.null(pts)) return(NULL)
  sb <- score_from_points(db, pts); so <- score_from_points(dat, pts)
  map_dfr(KS, function(k) { sa <- ss_k(sb, db$bbs, k); so_k <- ss_k(so, dat$bbs, k)
    tibble(rep = i, threshold = k, opt_sens = sa[["sens"]] - so_k[["sens"]], opt_spec = sa[["spec"]] - so_k[["spec"]]) })
}) %>% group_by(threshold) %>% summarise(opt_sens = mean(opt_sens), opt_spec = mean(opt_spec), .groups = "drop")

threshold_choice <- apparent_k %>% left_join(opt_k, by = "threshold") %>%
  mutate(sens_corr = sens_app - opt_sens, spec_corr = spec_app - opt_spec, youden_corr = sens_corr + spec_corr - 1) %>%
  group_by(score) %>% mutate(chosen = youden_corr == max(youden_corr)) %>% ungroup()
chosen_k <- threshold_choice %>% filter(chosen) %>% select(score, threshold, sens_corr, spec_corr, youden_corr)

# point stability (from 04 bootstrap)
point_stability <- map_dfr(ITEMS, function(it) { v <- boot[[paste0("pts_", it)]]
  tibble(item = it, legacy = c(SCORE_A_LEGACY, dental = 0)[it], rederived_V21 = pts_v21[it],
         pct_equal_rederived = 100 * mean(v == pts_v21[it]), pct_equal_legacy = 100 * mean(v == c(SCORE_A_LEGACY, dental = 0)[it]),
         pct_ge1 = 100 * mean(v >= 1)) })

# head-to-head summary of the two candidate scores at their chosen thresholds
k_v21 <- chosen_k$threshold[chosen_k$score == "rederived_V21"]; k_leg <- chosen_k$threshold[chosen_k$score == "legacy_scoreA"]
gA <- dat %>% filter(is.na(bbs_cat) | bbs_cat == "A") %>% mutate(bbs = as.integer(bbs_cat %in% "A"))
oth <- filter(dat, high_site == 0)
h2h <- map_dfr(names(SCORES), function(nm) { v <- SCORES[[nm]]; k <- if (nm == "rederived_V21") k_v21 else k_leg
  tibble(score = nm, max_points = if (nm == "rederived_V21") sum(pts_v21) else 13, threshold = k,
         auc_apparent = auc_fast(dat[[v]], dat$bbs),
         auc_corrected = auc_fast(dat[[v]], dat$bbs) - mean(boot$opt_auc_score),
         bbs_identified = sprintf("%d/38", sum(dat[[v]][dat$bbs == 1] >= k)),
         spec = 100 * mean(dat[[v]][dat$bbs == 0] < k),
         catA_identified = sprintf("%d/22", sum(gA[[v]][gA$bbs == 1] >= k)),
         auc_catA = auc_fast(gA[[v]], gA$bbs),
         auc_other_sites = auc_fast(oth[[v]], oth$bbs),
         other_sites_identified = sprintf("%d/7", sum(oth[[v]][oth$bbs == 1] >= k)),
         lt18_identified = sprintf("%d/23", sum(dat[[v]][dat$bbs == 1 & dat$age18 %in% "<18"] >= k)),
         ge18_identified = sprintf("%d/15", sum(dat[[v]][dat$bbs == 1 & dat$age18 %in% ">=18"] >= k)))
})

write_xlsx(list(threshold_choice = threshold_choice, chosen_threshold = chosen_k, point_stability = point_stability,
                head_to_head = h2h), file.path(OUT_DIR, "tables", "04c_score_selection.xlsx"))
message("04c score selection done: re-derived points ", paste(names(pts_v21), pts_v21, sep = "=", collapse = " "),
        "; chosen thresholds ", paste(chosen_k$score, chosen_k$threshold, collapse = ", "))
