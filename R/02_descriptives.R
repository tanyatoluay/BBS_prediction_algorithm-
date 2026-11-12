# =============================================================================
# 02 — Descriptive tables: Table 1, Table 2, Table 5, S5, S13, S16, adult ophthalmology
# =============================================================================
dat <- readRDS(file.path(OUT_DIR, "wegio_analysis_dataset.rds"))
set.seed(SEED_BOOT)

BIN_VARS <- c(age_lt18 = "Age <18 years", female = "Female sex", white = "White ethnicity",
  clin_dx = "Clinical diagnosis of BBS before enrolment", retinal = "Retinal dystrophy",
  visual = "Visual loss or optic atrophy", poly = "Polydactyly", syn = "Syndactyly", brachy = "Brachydactyly",
  renal = "Renal anomaly or dysfunction", hypog = "Hypogonadism", hypog_genital = "Hypogonadism or micropenis",
  cog = "Cognitive impairment", id_dis = "Intellectual disability", speech = "Speech delay or disorder",
  motor = "Delayed motor development", nd = "Neurodevelopmental impairment (any)", ataxia = "Ataxia",
  seizures = "Seizures", hyperphagia = "Hyperphagia", dental = "Dental crowding or hypodontia",
  hypert = "Arterial hypertension", cardiomyo = "Cardiomyopathy", hyperpig = "Hyperpigmented skin",
  dom3 = "≥3 affected organ-system domains", beales = "Beales criteria met", ern = "ERN criteria met")
dat$age_lt18 <- as.integer(dat$age18 == "<18")   # NA age -> NA (primary split)

exact_p <- function(x, g) {   # Fisher–Freeman–Halton exact test (r x 2 tables; computed exactly, no simulation)
  ok <- !is.na(x); tb <- table(g[ok, drop = TRUE], x[ok])
  if (ncol(tb) < 2) return(NA_real_)
  fisher.test(tb, workspace = 2e8)$p.value
}
kw_p <- function(x, g) { ok <- !is.na(x); kruskal.test(x[ok], droplevels(g[ok]))$p.value }
fmt_p <- function(p) ifelse(is.na(p), "—", ifelse(p < 0.001, "<0.001", ifelse(p < 0.01, sprintf("%.3f", p), sprintf("%.2f", p))))

g3 <- dat$grp4 %in% 1:3
tab1_bin <- map_dfr(names(BIN_VARS), function(v) {
  x <- dat[[v]]
  cells <- map_chr(levels(dat$group), ~ { xx <- x[dat$group == .x]; if (v == "hyperphagia") sprintf("%s [n=%d]", npct(sum(xx, na.rm = TRUE), sum(!is.na(xx))), sum(!is.na(xx))) else npct(sum(xx, na.rm = TRUE), length(xx)) })
  tibble(Characteristic = BIN_VARS[[v]], !!!setNames(as.list(cells), levels(dat$group)),
         Overall = if (v == "hyperphagia") sprintf("%s [n=%d]", npct(sum(x, na.rm = TRUE), sum(!is.na(x))), sum(!is.na(x))) else npct(sum(x, na.rm = TRUE), length(x)),
         p3 = exact_p(x[g3], dat$group[g3]), p4 = exact_p(x, dat$group))
})
cont <- list(
  "Age at enrolment, years"          = list(v = "age_y", f = function(g) TRUE, d = 1),
  "BMI-SDS (AGA), children"          = list(v = "bmi_sds", f = function(g) TRUE, d = 2),
  "BMI, kg/m², participants ≥16 years" = list(v = "bmi_calc", f = function(g) dat$age16 %in% ">=16", d = 1),
  "Recorded clinical features, n"    = list(v = "n_features", f = function(g) TRUE, d = 0),
  "Affected organ-system domains, n" = list(v = "n_domains", f = function(g) TRUE, d = 0))
tab1_cont <- imap_dfr(cont, function(s, lab) {
  sel <- rep_len(s$f(), nrow(dat)); x <- ifelse(sel, dat[[s$v]], NA)
  cells <- map_chr(levels(dat$group), ~ { xx <- x[dat$group == .x]; paste0(med_iqr(xx, s$d), " [n=", sum(!is.na(xx)), "]") })
  tibble(Characteristic = lab, !!!setNames(as.list(cells), levels(dat$group)),
         Overall = paste0(med_iqr(x, s$d), " [n=", sum(!is.na(x)), "]"),
         p3 = kw_p(x[g3], dat$group[g3]), p4 = kw_p(x, dat$group))
})
table1 <- bind_rows(tab1_cont, tab1_bin) %>% mutate(across(c(p3, p4), fmt_p))

# ---- Table 2: yields by age stratum ----------------------------------------------------------
yield_rows <- list(BBS = "bbs", `Non-BBS obesity-associated` = quote(grp4 == 2), `Non-obesity` = quote(grp4 == 3),
                   `Obesity-associated genetic disease` = "obesity_gd", `Genetic disease` = "genetic",
                   `No qualifying primary genetic diagnosis` = quote(grp4 == 4))
yv <- function(e) if (is.character(e)) dat[[e]] else as.integer(eval(e, dat))
table2 <- imap_dfr(yield_rows, function(e, lab) {
  y <- yv(e)
  strat <- function(a) { s <- !is.na(a); tb <- table(a[s], y[s])
    list(cells = map_chr(levels(a), ~ { yy <- y[a %in% .x]; sprintf("%d (%.1f; %s)", sum(yy), 100 * mean(yy),
            paste(sprintf("%.1f", wilson(sum(yy), length(yy))[2:3]), collapse = "–")) }),
         p = fisher.test(tb)$p.value) }
  s18 <- strat(dat$age18); s16 <- strat(dat$age16)
  tibble(Group = lab, All = sprintf("%d (%s)", sum(y), fmt_ci(sum(y), length(y))),
         `<18` = s18$cells[1], `>=18` = s18$cells[2], p18 = round(s18$p, 2),
         `<16` = s16$cells[1], `>=16` = s16$cells[2], p16 = round(s16$p, 2))
})

# ---- adult ophthalmology / protocol deviation -----------------------------------------------------
ad <- filter(dat, age18 == ">=18"); ch <- filter(dat, age18 == "<18")
ophth_summary <- tibble(
  metric = c("adults", "adults with ophthalmological feature", "  retinal dystrophy", "  visual loss only",
             "adults without ophthalmological feature", "  of whom protocol deviation", "children with ophth feature",
             "children with protocol deviation", "BBS in adults with ophth feature", "BBS in adults without",
             "Fisher p (BBS by ophth, adults)", "genetic disease % with ophth", "genetic disease % without"),
  value = c(nrow(ad), sum(ad$ophth), sum(ad$retinal), sum(ad$visual & !ad$retinal), sum(!ad$ophth),
            sum(ad$pd[ad$ophth == 0]), sum(ch$ophth), sum(ch$pd),
            sprintf("%d/%d (%.1f%%)", sum(ad$bbs[ad$ophth == 1]), sum(ad$ophth), 100 * mean(ad$bbs[ad$ophth == 1])),
            sprintf("%d/%d (%.1f%%)", sum(ad$bbs[ad$ophth == 0]), sum(!ad$ophth), 100 * mean(ad$bbs[ad$ophth == 0])),
            sprintf("%.4f", fisher.test(table(ad$ophth, ad$bbs))$p.value),
            sprintf("%.1f", 100 * mean(ad$genetic[ad$ophth == 1])), sprintf("%.1f", 100 * mean(ad$genetic[ad$ophth == 0]))))

# ---- Table 5: within-group age comparison (<16 vs >=16) --------------------------------------------
T5_VARS <- c("retinal", "visual", "poly", "synbrachy", "renal", "hypog_genital", "nd", "id_dis", "hyperphagia",
             "dental", "hypert", "dom3", "beales", "clin_dx")
age_compare <- function(split) map_dfr(levels(dat$group), function(gname) map_dfr(T5_VARS, function(v) {
  s <- dat %>% filter(group == gname, !is.na(.data[[split]]))
  a <- na.omit(s[[v]][s[[split]] == levels(s[[split]])[1]]); b <- na.omit(s[[v]][s[[split]] == levels(s[[split]])[2]])
  tibble(group = gname, feature = v, young = sprintf("%d/%d", sum(a), length(a)), old = sprintf("%d/%d", sum(b), length(b)),
         p = fisher.test(matrix(c(sum(a), length(a) - sum(a), sum(b), length(b) - sum(b)), 2))$p.value)
}))
table5_16 <- age_compare("age16"); table5_18 <- age_compare("age18")

bbs_noret_age <- with(filter(dat, bbs == 1), wilcox.test(age_y ~ retinal, exact = FALSE))
bbs_age_by_ret <- dat %>% filter(bbs == 1) %>% group_by(retinal) %>% summarise(median_age = median(age_y))

# ---- Table S5: pairwise Fisher + BH --------------------------------------------------------------------
S5_VARS <- c("female", "retinal", "visual", "poly", "syn", "brachy", "renal", "hypog", "hypog_genital", "cog",
             "id_dis", "speech", "motor", "nd", "ataxia", "seizures", "hyperphagia", "dental", "hypert",
             "cardiomyo", "hyperpig", "dom3", "beales", "ern")          # 24 features
pairs <- list(c(1, 2), c(1, 3), c(2, 3), c(1, 4), c(2, 4), c(3, 4))
tableS5 <- map_dfr(pairs, function(pr) {
  s <- filter(dat, grp4 %in% pr)
  map_dfr(S5_VARS, function(v) {
    tb <- table(factor(s$grp4, levels = pr), factor(s[[v]], levels = c(1, 0)))
    ft <- fisher.test(tb)
    a <- tb[1, 1]; b <- tb[1, 2]; c <- tb[2, 1]; d_ <- tb[2, 2]
    tibble(comparison = paste(levels(dat$group)[pr], collapse = " vs "), feature = v,
           first = sprintf("%d/%d", a, a + b), second = sprintf("%d/%d", c, c + d_),
           OR = (a + .5) * (d_ + .5) / ((b + .5) * (c + .5)),       # Haldane–Anscombe corrected odds ratio
           OR_cond_mle = unname(ft$estimate), p = ft$p.value)
  }) %>% mutate(p_adj_within = p.adjust(p, "BH"))
}) %>% mutate(p_adj_all144 = p.adjust(p, "BH"))

# ---- Table S13: three-state documentation ---------------------------------------------------------------
tableS13 <- map_dfr(c("retinal", "renal", "hypog_genital", "poly", "nd"), function(f)
  dat %>% filter(!is.na(age18)) %>% count(feature = f, group, age18, state = .data[[paste0(f, "_3")]]) %>%
    pivot_wider(names_from = state, values_from = n, values_fill = 0))

# ---- Table S16: sites -----------------------------------------------------------------------------------------
tableS16 <- dat %>% group_by(site) %>% summarise(n = n(), BBS = sum(bbs), nonBBS_obesity = sum(grp4 == 2),
                                                  non_obesity = sum(grp4 == 3), genetic = sum(genetic)) %>% arrange(desc(BBS))
site_summary <- tibble(
  n_sites = n_distinct(dat$site), sites_with_BBS = sum(tableS16$BBS > 0),
  high_sites_n = sum(dat$high_site), high_sites_BBS = sum(dat$bbs[dat$high_site == 1]),
  other_n = sum(!dat$high_site), other_BBS = sum(dat$bbs[dat$high_site == 0]),
  other_genetic = sum(dat$genetic[dat$high_site == 0]))

write_xlsx(list(table1 = table1, table2 = table2, ophth_adults = ophth_summary, table5_age16 = table5_16,
                table5_age18 = table5_18, tableS5 = tableS5, tableS13 = tableS13, tableS16 = tableS16,
                site_summary = site_summary),
           file.path(OUT_DIR, "tables", "02_descriptives.xlsx"))
message("02 descriptives done")
