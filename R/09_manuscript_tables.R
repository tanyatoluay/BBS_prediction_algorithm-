# =============================================================================
# 09 — Every data-derived manuscript table, rebuilt in the manuscript layout (character cells,
#      mid-dot decimals). Runs after 01–08 in the same session. Output: output/manuscript_tables.xlsx
#      Tables whose content is a genetic/biological classification (table 3, S3, S10, S2b, S1 mapping,
#      genotype column of S4) are inputs, not statistics, and are not rebuilt here.
# =============================================================================
md  <- function(x, d = 1) ifelse(is.na(x), "—", chartr(".", "·", formatC(x, format = "f", digits = d)))
pct <- function(k, n, d = 1) md(100 * k / n, d)
npc <- function(k, n, d = 1) sprintf("%d (%s%%)", k, pct(k, n, d))
wci <- function(k, n) { w <- wilson(k, n); sprintf("%s (%s–%s)", md(w[1]), md(w[2]), md(w[3])) }
p2  <- function(p) ifelse(p < 0.001, "<0·001", ifelse(p < 0.01, md(p, 3), md(p, 2)))           # tables 1, 5
p4  <- function(p) ifelse(p < 0.0001, "<0·0001", ifelse(p < 0.001, md(p, 4), ifelse(p < 0.01, md(p, 3), md(p, 2))))  # tables 6, S5, S6
q7  <- function(v, d) { q <- quantile(v, c(.5, .25, .75), na.rm = TRUE, type = 7)
  f <- function(z) if (d == 0 && z != round(z)) md(z, 1) else md(z, d); sprintf("%s (%s–%s)", f(q[1]), f(q[2]), f(q[3])) }
d <- dat; G <- levels(d$group); MT <- list()
LAB <- c(retinal = "Retinal dystrophy", visual = "Visual loss or optic atrophy", poly = "Polydactyly", syn = "Syndactyly",
         brachy = "Brachydactyly", synbrachy = "Syndactyly or brachydactyly", renal = "Renal anomaly or dysfunction", hypog = "Hypogonadism",
         hypog_genital = "Hypogonadism or micropenis", cog = "Cognitive impairment", id_dis = "Intellectual disability",
         speech = "Speech delay or disorder", motor = "Delayed motor development", nd = "Neurodevelopmental impairment", ataxia = "Ataxia",
         seizures = "Seizures", hyperphagia = "Hyperphagia", dental = "Dental crowding or hypodontia", hypert = "Arterial hypertension",
         cardiomyo = "Cardiomyopathy", hyperpig = "Hyperpigmented skin", dom3 = "≥3 organ-system domains", beales = "Beales criteria met",
         ern = "ERN criteria met", female = "Female sex", male = "Male sex", clin_dx = "Clinical BBS diagnosis before enrolment")
mat <- function(rows) do.call(rbind, lapply(rows, function(r) as.character(r)))

# ---- Table 1 ---------------------------------------------------------------------------------------------
g3 <- d$grp4 < 4
t1_bin <- function(lab, x, assessed_only = FALSE) {
  cell <- function(xx, n_all) if (assessed_only) sprintf("%s [n=%d]", npc(sum(xx, na.rm = TRUE), sum(!is.na(xx))), sum(!is.na(xx))) else npc(sum(xx, na.rm = TRUE), n_all)
  cells <- sapply(G, function(g) cell(x[d$group == g], sum(d$group == g)))
  c(lab, cells, cell(x, nrow(d)), p2(exact_p(x[g3], d$group[g3])), p2(exact_p(x, d$group))) }
t1_con <- function(lab, x, dd, shown_n = TRUE) { cells <- sapply(G, function(g) { v <- x[d$group == g]
  paste0(q7(v, dd), if (shown_n) sprintf(" [n=%d]", sum(!is.na(v))) else "") })
  c(lab, cells, paste0(q7(x, dd), if (shown_n) sprintf(" [n=%d]", sum(!is.na(x))) else ""),
    p2(kw_p(x[g3], d$group[g3])), p2(kw_p(x, d$group))) }
bmi16 <- ifelse(d$age16 %in% ">=16", d$bmi_calc, NA)
MT$table1 <- mat(list(
  t1_con("Age at enrolment, years, median (IQR)", d$age_y, 1, FALSE),
  t1_bin("Age <18 years", as.integer(d$age18 == "<18")), t1_bin("Female sex", d$female), t1_bin("White ethnicity", d$white),
  t1_con("BMI-SDS (AGA reference), children", d$bmi_sds, 2), t1_con("BMI, kg/m², participants ≥16 years", bmi16, 1),
  t1_bin("Clinical diagnosis of BBS recorded before enrolment", d$clin_dx),
  t1_bin("Retinal dystrophy (rod–cone dystrophy, retinal degeneration, retinitis pigmentosa)", d$retinal),
  t1_bin("Visual loss or optic atrophy", d$visual), t1_bin("Polydactyly", d$poly), t1_bin("Syndactyly", d$syn),
  t1_bin("Brachydactyly", d$brachy), t1_bin("Renal anomaly or dysfunction", d$renal), t1_bin("Hypogonadism", d$hypog),
  t1_bin("Hypogonadism or micropenis", d$hypog_genital), t1_bin("Cognitive impairment", d$cog), t1_bin("Intellectual disability", d$id_dis),
  t1_bin("Speech delay or disorder", d$speech), t1_bin("Delayed motor development", d$motor),
  t1_bin("Neurodevelopmental impairment (any)", d$nd), t1_bin("Ataxia", d$ataxia), t1_bin("Seizures", d$seizures),
  t1_bin("Hyperphagia", d$hyperphagia, assessed_only = TRUE), t1_bin("Dental crowding or hypodontia", d$dental), t1_bin("Arterial hypertension", d$hypert),
  t1_bin("Cardiomyopathy", d$cardiomyo), t1_bin("Hyperpigmented skin", d$hyperpig),
  t1_con("Recorded clinical features, n, median (IQR)", d$n_features, 0, FALSE),
  t1_con("Affected organ-system domains, n, median (IQR)", d$n_domains, 0, FALSE),
  t1_bin("≥3 affected organ-system domains", d$dom3), t1_bin("Beales criteria met", d$beales), t1_bin("ERN criteria met", d$ern)))

# ---- Table 2 (upper) ---------------------------------------------------------------------------------------------
y2 <- list("BBS" = d$bbs, "Non-BBS obesity-associated" = as.integer(d$grp4 == 2), "Non-obesity" = as.integer(d$grp4 == 3),
           "Obesity-associated genetic disease (BBS + non-BBS)" = d$obesity_gd, "Genetic disease (all three groups)" = d$genetic,
           "No qualifying primary genetic diagnosis" = as.integer(d$grp4 == 4))
cell_ci <- function(k, n) { w <- wilson(k, n); sprintf("%d (%s; %s–%s)", k, md(w[1]), md(w[2]), md(w[3])) }
fpp <- function(y, a) { s <- !is.na(a); md(fisher.test(table(a[s], y[s]))$p.value, 2) }
MT$table2 <- mat(lapply(names(y2), function(lab) { y <- y2[[lab]]
  a <- d$age18; b <- d$age16
  c(lab, cell_ci(sum(y), length(y)),
    sprintf("%d (%s)", sum(y[a %in% "<18"]), pct(sum(y[a %in% "<18"]), sum(a %in% "<18"))),
    sprintf("%d (%s)", sum(y[a %in% ">=18"]), pct(sum(y[a %in% ">=18"]), sum(a %in% ">=18"))), fpp(y, a),
    cell_ci(sum(y[b %in% "<16"]), sum(b %in% "<16")), cell_ci(sum(y[b %in% ">=16"]), sum(b %in% ">=16")), fpp(y, b)) }))

# ---- Table 4 -------------------------------------------------------------------------------------------------------
crit <- c("Clinical diagnosis before enrolment" = "clin_dx", "ERN criteria" = "ern", "Beales criteria" = "beales",
          setNames("scoreA3", paste0("Score A ≥", THR_LO)), setNames("scoreAhi", paste0("Score A ≥", THR_HI)))
strata <- list("All" = rep(TRUE, nrow(d)), "<18 years" = d$age18 %in% "<18", "≥18 years" = d$age18 %in% ">=18",
               "<16 years" = d$age16 %in% "<16", "≥16 years" = d$age16 %in% ">=16")
MT$table4 <- mat(unlist(lapply(names(crit), function(cl) lapply(names(strata), function(sl) {
  s <- strata[[sl]]; x <- d[[crit[[cl]]]][s]; y <- d$bbs[s]
  tp <- sum(x & y); fn <- sum(!x & y); fp_ <- sum(x & !y); tn <- sum(!x & !y)
  c(if (sl == "All") cl else "", sl, sprintf("%d/%d", tp, tp + fn), wci(tp, tp + fn), wci(tn, tn + fp_), tp + fp_,
    md(100 * tp / (tp + fp_)), md(100 * tn / (tn + fn))) })), recursive = FALSE))

# ---- Table 5 (<18 vs >=18) ---------------------------------------------------------------------------------------------
T5V <- c("retinal", "visual", "poly", "synbrachy", "renal", "hypog_genital", "nd", "id_dis", "hyperphagia", "dental", "hypert", "dom3", "beales", "clin_dx")
np0 <- function(k, n) sprintf("%d (%s%%)", k, md(100 * k / n, 0))
MT$table5 <- mat(lapply(T5V, function(v) c(LAB[[v]], unlist(lapply(G, function(g) { s <- filter(d, group == g, !is.na(age18))
  a <- na.omit(s[[v]][s$age18 == "<18"]); b <- na.omit(s[[v]][s$age18 == ">=18"])
  fmt <- if (v == "hyperphagia") function(k, n) sprintf("%d/%d (%s%%)", k, n, md(100 * k / n, 0)) else np0
  c(fmt(sum(a), length(a)), fmt(sum(b), length(b)), p2(fisher.test(matrix(c(sum(a), length(a) - sum(a), sum(b), length(b) - sum(b)), 2))$p.value)) })))))

# ---- Table 6 upper and lower -------------------------------------------------------------------------------------------
t6 <- table6_model %>% filter(term != "(Intercept)")
t6lab <- c(poly = "Polydactyly", retinal = "Retinal dystrophy", hypog_genital = "Hypogonadism or genital hypoplasia",
           renal = "Renal anomaly or dysfunction", synbrachy = "Syndactyly or brachydactyly", nd = "Neurodevelopmental impairment",
           dental = "Dental anomaly", age10 = "Age (per 10 years)", male = "Male sex")
MT$table6_upper <- mat(lapply(seq_len(nrow(t6)), function(i) { r <- t6[i, ]; it <- r$term
  c(t6lab[[it]], md(r$beta, 2), sprintf("%s (%s–%s)", md(r$OR, 1), md(r$OR_pl_lo, 1), md(r$OR_pl_hi, 1)), p4(r$p_pl),
    if (it %in% names(SCORE_A_POINTS)) SCORE_A_POINTS[[it]] else if (it == "dental") "not scored" else "—",
    if (it %in% names(SCORE_A_LEGACY)) SCORE_A_LEGACY[[it]] else if (it == "dental") "0" else "—") }))
MT$table6_upper[, 2] <- sub("^-", "-", MT$table6_upper[, 2])
MT$table6_lower <- mat(lapply(1:sum(SCORE_A_POINTS), function(k) { x <- d$scoreA >= k; y <- d$bbs
  tp <- sum(x & y); fn <- sum(!x & y); fp_ <- sum(x & !y); tn <- sum(!x & !y)
  c(paste0("≥", k), tp + fp_, pct(tp + fp_, nrow(d)), tp, wci(tp, tp + fn), wci(tn, tn + fp_), wci(tp, tp + fp_), md(100 * tn / (tn + fn))) }))

# ---- Table 7 (bands) ----------------------------------------------------------------------------------------------------
MT$table7 <- mat(lapply(seq_len(nrow(table7)), function(i) { r <- table7[i, ]
  c(as.character(r$band), r$n, r$bbs, wci(r$bbs, r$n), r$nonBBS_ob, wci(r$bbs + r$nonBBS_ob, r$n), r$nonob,
    wci(r$bbs + r$nonBBS_ob + r$nonob, r$n), r$none) }))

# ---- Table 8 ------------------------------------------------------------------------------------------------------------
gA <- d %>% filter(is.na(bbs_cat) | bbs_cat == "A") %>% mutate(bbs = as.integer(bbs_cat %in% "A"))
pops <- list("All study-assigned BBS (n=38)" = d, "Molecularly confirmed: biallelic P/LP (category A); other BBS excluded" = gA,
  "Molecular rule met (A + B1); below-rule BBS excluded" = filter(d, is.na(bbs_cat) | bbs_cat %in% c("A", "B1")),
  "Molecularly confirmed (A) vs all others" = mutate(d, bbs = as.integer(bbs_cat %in% "A")),
  "Index participants only (n=983)" = filter(d, index == 1), "Excluding the two sites with most BBS cases" = filter(d, high_site == 0),
  "Participants <18 years (paediatric eligibility)" = filter(d, age18 %in% "<18"), "Participants ≥18 years (adult eligibility)" = filter(d, age18 %in% ">=18"))
sv <- function(p, v, spec = TRUE) { tp <- sum(p[[v]] & p$bbs == 1); s <- sprintf("%d (%s)", tp, pct(tp, sum(p$bbs)))
  if (spec) paste0(s, "; ", md(100 * mean(p[[v]][p$bbs == 0] == 0))) else s }
MT$table8 <- mat(lapply(names(pops), function(nm) { p <- pops[[nm]]; a <- auc_boot_ci(p$scoreA, p$bbs)
  c(nm, nrow(p), sum(p$bbs), sprintf("%s (%s–%s)", md(a[1], 3), md(a[2], 3), md(a[3], 3)), sv(p, "scoreA3", FALSE), sv(p, "scoreAhi"),
    sv(p, "beales"), sv(p, "ern")) }))

# ---- Table S2 -----------------------------------------------------------------------------------------------------------
bb <- filter(d, bbs == 1); nb <- filter(d, bbs == 0)
s2c <- function(v) c(sum(d[[v]]), sprintf("%d (%s%%)", sum(bb[[v]]), pct(sum(bb[[v]]), 38)), sum(bb[[v]][bb$age16 %in% "<16"]),
                     sum(bb[[v]][bb$age16 %in% ">=16"]), sum(nb[[v]]), sum(nb[[v]][nb$grp4 == 4]), md(100 * mean(nb[[v]] == 0)))
a_ <- s2c("beales_export"); b_ <- s2c("beales")
MT$tableS2 <- cbind(c("Beales criteria met, all participants", "Beales criteria met, BBS (n=38)", "   BBS <16 years (n=23)", "   BBS ≥16 years (n=15)",
  "Beales criteria met without BBS", "   of which no genetic disease", "Specificity, %", "ERN criteria met, all participants", "ERN criteria met, BBS"),
  c(a_, sum(d$ern_export), sum(bb$ern_export)), c(b_, sum(d$ern), sum(bb$ern)))

# ---- Table S4 (all columns except genotype; matched by participant) --------------------------------------------------------
codes <- function(r) paste(c(if (r$poly) "Po", if (r$retinal) "Re", if (r$renal) "Rn", if (r$hypog_genital) "Ge", if (r$synbrachy) "SB",
                             if (r$nd) "ND", if (r$dental) "De", if (isTRUE(r$hyperphagia_any == 1)) "Hp"), collapse = " ")
MT$tableS4 <- mat(lapply(seq_len(nrow(bb)), function(i) { r <- bb[i, ]
  c(r$id, md(r$age_y, 1), ifelse(r$female == 1, "F", "M"), NA, r$bbs_cat, ifelse(r$beales == 1, "yes", "no"), ifelse(r$ern == 1, "yes", "no"),
    r$scoreA, ifelse(r$clin_dx == 1, "yes", "no"), codes(r)) }))

# ---- Table S5 ----------------------------------------------------------------------------------------------------------------
cmp <- c("BBS vs Non-BBS obesity-associated" = "BBS vs Non-BBS obesity", "BBS vs Non-obesity" = "BBS vs Non-obesity",
         "Non-BBS obesity-associated vs Non-obesity" = "Non-BBS obesity vs Non-obesity", "BBS vs No qualifying diagnosis" = "BBS vs No genetic disease",
         "Non-BBS obesity-associated vs No qualifying diagnosis" = "Non-BBS obesity vs No genetic disease",
         "Non-obesity vs No qualifying diagnosis" = "Non-obesity vs No genetic disease")
s5 <- tableS5 %>% filter(p_adj_within < 0.10)
MT$tableS5 <- mat(lapply(seq_len(nrow(s5)), function(i) { r <- s5[i, ]
  c(cmp[[r$comparison]], LAB[[r$feature]], r$first, r$second, md(r$OR, 2), p4(r$p), p4(r$p_adj_within)) }))

# ---- Table S6 upper / lower ------------------------------------------------------------------------------------------------
F16 <- c("male", "retinal", "visual", "poly", "synbrachy", "renal", "hypog_genital", "cog", "id_dis", "speech", "motor", "ataxia", "seizures",
         "hyperphagia", "dental", "hypert")
s6 <- function(case_grp, lab) { s <- filter(d, grp4 %in% c(case_grp, 4)); y <- as.integer(s$grp4 %in% case_grp)
  r <- map_dfr(F16, function(f) { x <- s[[f]]; ok <- !is.na(x); x <- x[ok]; yy <- y[ok]
    a <- sum(x & yy); b <- sum(x & !yy); c <- sum(!x & yy); dd <- sum(!x & !yy)
    tibble(f = f, n = a + b, a = a, ppv = md(100 * a / (a + b)), sens = md(100 * a / sum(yy)),
           OR = md((a + .5) * (dd + .5) / ((b + .5) * (c + .5)), 2), p = fisher.test(matrix(c(a, b, c, dd), 2))$p.value) })
  r$padj <- p.adjust(r$p, "BH"); r <- arrange(r, p) %>% head(8)
  lapply(seq_len(nrow(r)), function(i) c(lab, LAB[[r$f[i]]], r$n[i], r$a[i], r$ppv[i], r$sens[i], r$OR[i], p4(r$padj[i]))) }
MT$tableS6_upper <- mat(c(s6(2, "Non-BBS obesity-associated vs no genetic disease"), s6(2:3, "Non-BBS genetic disease vs no genetic disease")))
cnt <- nb$nd + nb$hypog_genital + nb$id_dis
MT$tableS6_lower <- mat(lapply(0:3, function(k) { s <- cnt >= k
  c(paste0("≥", k), sum(s), pct(sum(s), nrow(nb)), sprintf("%d (%s%%)", sum(nb$genetic[s]), pct(sum(nb$genetic[s]), sum(nb$genetic))),
    pct(sum(nb$genetic[s]), sum(s)), sprintf("%d (%s%%)", sum(nb$grp4[s] == 2), pct(sum(nb$grp4[s] == 2), sum(nb$grp4 == 2))),
    pct(sum(nb$grp4[s] == 2), sum(s))) }))

# ---- Tables S7, S8, S11, S13, S14, S15, S16 -----------------------------------------------------------------------------------
s7r <- function(lab, p, sc = "scoreA") { x <- p[[sc]] >= THR_LO; y <- p$bbs; tp <- sum(x & y); fn <- sum(!x & y); fp_ <- sum(x & !y); tn <- sum(!x & !y)
  c(lab, nrow(p), sum(y), md(auc_fast(p[[sc]], y), 3), pct(tp, tp + fn), pct(tn, tn + fp_), pct(tp, tp + fp_), md(100 * tn / (tn + fn))) }
MT$tableS7 <- mat(list(
  s7r("Primary analysis (BBS group, n=38)", d),
  s7r("BBS meeting the molecular rule; 11 below-rule participants excluded", filter(d, !(bbs_cat %in% c("B2", "C")))),
  s7r("BBS meeting the molecular rule; 11 below-rule participants as non-cases", mutate(d, bbs = as.integer(bbs_cat %in% c("A", "B1")))),
  s7r("Participants without a clinical BBS diagnosis before enrolment", filter(d, clin_dx == 0)),
  s7r("Participants without recorded protocol deviation", filter(d, pd == 0)),
  s7r("Participants <16 years", filter(d, age16 %in% "<16")),
  s7r("Participants <16 years, retinal item removed", filter(d, age16 %in% "<16"), "scoreA_noret"),
  s7r("Participants ≥16 years", filter(d, age16 %in% ">=16")),
  s7r("Retinal item removed, all participants", d, "scoreA_noret"),
  s7r("Participants <6 years", filter(d, age_y < 6)),
  s7r("Not-assessed features counted as present", d, "scoreA_napres"),
  s7r("Participants with all score items assessed (complete case)", filter(d, scoreA_complete))))
MT$tableS8 <- mat(unlist(lapply(c("<16", ">=16"), function(a) { p <- filter(d, age16 %in% a); lapply(c(2:6, 8), function(k) {
  x <- p$scoreA >= k; y <- p$bbs; tp <- sum(x & y); fn <- sum(!x & y); fp_ <- sum(x & !y); tn <- sum(!x & !y)
  c(ifelse(a == "<16", "<16 years", "≥16 years"), paste0("≥", k), tp + fp_, tp, pct(tp, tp + fn), pct(tn, tn + fp_), pct(tp, tp + fp_), md(100 * tn / (tn + fn))) }) }),
  recursive = FALSE))
catlab <- c(A = "A — biallelic P/LP (homozygous or compound heterozygous)", B1 = "B1 — P/LP and VUS-high second allele",
            B2 = "B2 — below molecular rule; Beales or ERN met", C = "C — below molecular rule; neither criterion met")
MT$tableS11 <- mat(lapply(names(catlab), function(k) { s <- filter(bb, bbs_cat == k)
  c(catlab[[k]], nrow(s), sum(s$age18 == "<18"), sum(s$beales), sum(s$ern), sum(s$scoreAhi), sum(s$clin_dx), NA) }))
s13lab <- c(retinal = "Retinal dystrophy", renal = "Renal anomaly or dysfunction", hypog_genital = "Hypogonadism or micropenis",
            poly = "Polydactyly", nd = "Neurodevelopmental impairment")
glab <- c("BBS", "Non-BBS obesity", "Non-obesity", "No qualifying diagnosis")
MT$tableS13 <- mat(unlist(lapply(names(s13lab), function(f) unlist(lapply(1:4, function(gi) lapply(c("<18", ">=18"), function(a) {
  s <- d[d$grp4 == gi & d$age18 %in% a, ]; st <- s[[paste0(f, "_3")]]
  c(s13lab[[f]], glab[gi], ifelse(a == "<18", "<18", "≥18"), nrow(s), sum(st == "present"), sum(st == "absent"), sum(st == "not_assessed")) })),
  recursive = FALSE)), recursive = FALSE))
s14lab <- t6lab[ITEMS]
MT$tableS14 <- mat(lapply(seq_len(nrow(tableS14)), function(i) { r <- tableS14[i, ]
  c(s14lab[[r$item]], r$scoreA, md(r$pct_ge1), sprintf("%d (%d–%d)", as.integer(r$median), as.integer(r$p2.5), as.integer(r$p97.5)), md(r$pct_equal_scoreA)) }))
MT$tableS15 <- mat(lapply(seq_len(nrow(tableS15)), function(i) c(md(100 * tableS15$prev[i]), md(tableS15$ppv[i]), md(tableS15$npv[i]))))
st <- d %>% group_by(site) %>% summarise(n = n(), b = sum(bbs), o = sum(grp4 == 2), no = sum(grp4 == 3)) %>% filter(site %in% c("10", "39", "16"))
st <- st[match(c("10", "39", "16"), st$site), ]
MT$tableS16 <- rbind(mat(lapply(seq_len(nrow(st)), function(i) c(st$site[i], st$n[i], st$b[i], st$o[i], st$no[i]))),
                     c("All sites except 10 and 39", sum(d$high_site == 0), sum(d$bbs[d$high_site == 0]), "—", "—"))

write_xlsx(lapply(MT, function(m) as.data.frame(m, stringsAsFactors = FALSE)), file.path(OUT_DIR, "manuscript_tables.xlsx"))
saveRDS(MT, file.path(OUT_DIR, "manuscript_tables.rds"))
message("09 manuscript tables written: ", length(MT), " tables")
