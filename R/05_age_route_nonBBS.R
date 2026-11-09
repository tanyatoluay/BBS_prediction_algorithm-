# =============================================================================
# 05 — Age strata by documented eligibility route; non-BBS prediction (table S6)
# =============================================================================
dat <- readRDS(file.path(OUT_DIR, "wegio_analysis_dataset.rds"))

# ---- 1. Adults by documented eCRF route ------------------------------------------------------------
# "adult section only"  = enrolled via adult criteria (obesity + rod–cone dystrophy / visual loss)
# "both sections"       = paediatric checklist (obesity <6 y + feature) completed in an adult -> protocol deviation
route_tab <- dat %>% filter(!is.na(age18)) %>% group_by(age18, route) %>%
  summarise(n = n(), protocol_dev = sum(pd), early_obesity_lt6 = sum(ad_obese_lt6 %in% "yes" | !is.na(ck_poly)),
            spec_retinal_nonBBS = round(100 * mean(retinal[bbs == 0] == 0), 1),
            BBS_pct = round(100 * mean(bbs), 1), genetic_pct = round(100 * mean(genetic), 1),
            BBS_site10_39 = sum(bbs & high_site),
            n_ophth = sum(ophth), n_retinal = sum(retinal), n_BBS = sum(bbs), n_genetic = sum(genetic), .groups = "drop")

ad <- filter(dat, age18 == ">=18", route %in% c("adult section only", "both sections"))
route_tests <- tibble(
  comparison = c("BBS: adult-route vs paediatric-route adults", "Genetic disease: adult-route vs paediatric-route adults"),
  p_fisher = c(fisher.test(table(ad$route, ad$bbs))$p.value, fisher.test(table(ad$route, ad$genetic))$p.value))

# Score A within routes (retinal item is an entry criterion in the adult route)
route_auc <- dat %>% filter(age18 == ">=18") %>% group_by(route) %>% filter(sum(bbs) > 0, sum(bbs) < n()) %>%
  summarise(auc_scoreA = auc_fast(scoreA, bbs), auc_noret = auc_fast(scoreA_noret, bbs),
            sens5 = sum(scoreA5 & bbs), n = n(), n_bbs = sum(bbs), .groups = "drop")

# ---- 2. Age comparison restricted to participants in whom the feature was assessed -------------------
# Hyperphagia exists only in the paediatric checklist: absent checklist = structurally missing, not "absent".
age_assessed <- map_dfr(c("retinal", "hyperphagia", "nd", "id_dis", "renal", "poly"), function(f) {
  d <- filter(dat, bbs == 1)
  assessed <- if (f == "hyperphagia") d$hyperphagia_assessed else d[[paste0(f, "_3")]] != "not_assessed"
  map_dfr(c("primary (not assessed = absent)", "assessed only"), function(m) {
    s <- if (m == "assessed only") d[assessed, ] else d
    y <- s[[f]][s$age18 == "<18"]; o <- s[[f]][s$age18 == ">=18"]
    tibble(feature = f, handling = m, lt18 = sprintf("%d/%d", sum(y), length(y)), ge18 = sprintf("%d/%d", sum(o), length(o)),
           p = fisher.test(matrix(c(sum(y), length(y) - sum(y), sum(o), length(o) - sum(o)), 2))$p.value)
  })
})

# ---- 3. Non-BBS genetic disease (table S6) ---------------------------------------------------------------------------
# 16 candidate features (assumed list; V25 does not enumerate them -> please confirm)
S6_FEATS <- c("male", "retinal", "visual", "poly", "synbrachy", "renal", "hypog_genital", "cog", "id_dis", "speech",
              "motor", "ataxia", "seizures", "hyperphagia", "dental", "hypert")
nb <- filter(dat, bbs == 0)
s6_fit <- function(y) {
  f <- reformulate(S6_FEATS, response = y)
  fit <- logistf(f, data = nb, pl = FALSE)
  lp0 <- drop(model.matrix(f, nb) %*% coef(fit)); app <- auc_fast(lp0, nb[[y]])
  set.seed(SEED_BOOT)
  opt <- map_dbl(1:300, function(i) { s <- sample(nrow(nb), replace = TRUE); db <- nb[s, ]
    fb <- logistf(f, data = db, pl = FALSE)
    auc_fast(drop(model.matrix(f, db) %*% coef(fb)), db[[y]]) - auc_fast(drop(model.matrix(f, nb) %*% coef(fb)), nb[[y]]) })
  tibble(outcome = y, cases = sum(nb[[y]]), auc_apparent = app, auc_corrected = app - mean(opt))
}
nb$nonbbs_gd <- nb$genetic; nb$nonbbs_ob <- as.integer(nb$grp4 == 2)
tableS6_model <- bind_rows(s6_fit("nonbbs_gd"), s6_fit("nonbbs_ob"))
cnt <- nb$nd + nb$hypog_genital + nb$id_dis
tableS6_count <- map_dfr(0:3, ~ tibble(count_ge = .x, n = sum(cnt >= .x), nonBBS_gd = sum(nb$genetic[cnt >= .x]),
  ppv_gd = round(100 * mean(nb$genetic[cnt >= .x]), 1), nonBBS_ob = sum(nb$grp4[cnt >= .x] == 2),
  ppv_ob = round(100 * mean(nb$grp4[cnt >= .x] == 2), 1)))

write_xlsx(list(route_by_age = route_tab, route_tests = route_tests, route_auc = route_auc,
                bbs_age_assessed_only = age_assessed, tableS6_model = tableS6_model, tableS6_count = tableS6_count),
           file.path(OUT_DIR, "tables", "05_age_route_nonBBS.xlsx"))
message("05 age/route/non-BBS done")
