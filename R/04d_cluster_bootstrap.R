# =============================================================================
# 04d — Site-cluster bootstrap
#   31 of 38 BBS cases come from two sites; participants within a site are not independent.
#   Whole sites are resampled with replacement (37 sites, 2000 resamples), and percentile 95% CIs
#   are computed for the AUC of Score A, sensitivity/specificity at THR_HI and THR_LO, and the
#   AUC difference vs the Beales feature count. Compared with the participant-level bootstrap.
# =============================================================================
N_CLUSTER_BOOT <- 2000
sites <- sort(unique(dat$site)); by_site <- split(seq_len(nrow(dat)), dat$site)

stat_fun <- function(idx) {
  y <- dat$bbs[idx]; s <- dat$scoreA[idx]
  if (sum(y) == 0 || sum(y) == length(y)) return(rep(NA_real_, 6))
  c(auc = auc_fast(s, y), diff = auc_fast(s, y) - auc_fast(dat$beales_count[idx], y),
    sens_hi = mean(s[y == 1] >= THR_HI), spec_hi = mean(s[y == 0] < THR_HI),
    sens_lo = mean(s[y == 1] >= THR_LO), spec_lo = mean(s[y == 0] < THR_LO))
}
set.seed(SEED_BOOT)
cb <- t(replicate(N_CLUSTER_BOOT, {
  pick <- sample(sites, length(sites), replace = TRUE)
  stat_fun(unlist(by_site[pick], use.names = FALSE))
}))
cb <- cb[complete.cases(cb), , drop = FALSE]

# participant-level (stratified) bootstrap for comparison
set.seed(SEED_BOOT); i1 <- which(dat$bbs == 1); i0 <- which(dat$bbs == 0)
pb <- t(replicate(N_CLUSTER_BOOT, stat_fun(c(rs(i1), rs(i0)))))

est <- stat_fun(seq_len(nrow(dat)))
cluster_ci <- tibble(
  statistic = names(est), estimate = est,
  participant_lo = apply(pb, 2, quantile, .025), participant_hi = apply(pb, 2, quantile, .975),
  cluster_lo = apply(cb, 2, quantile, .025), cluster_hi = apply(cb, 2, quantile, .975),
  n_cluster_resamples = nrow(cb))
cases_per_resample <- tibble(statistic = "BBS cases per cluster resample",
  median = median(replicate(500, sum(dat$bbs[unlist(by_site[sample(sites, length(sites), TRUE)])]))))

write_xlsx(list(cluster_vs_participant_CI = cluster_ci), file.path(OUT_DIR, "tables", "04d_cluster_bootstrap.xlsx"))
message("04d site-cluster bootstrap done (", nrow(cb), " valid resamples)")
