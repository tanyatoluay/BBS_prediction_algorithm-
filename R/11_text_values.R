# =============================================================================
# 11 — Every number in the running text (Summary, Research in context, Methods, Results, Discussion)
#      checked against R. Each claim is a text fragment generated from the data; the check looks for it
#      in the manuscript (accepted view). Numbers from genetic/classification inputs are listed as inputs.
#      A coverage step lists every number in those sections that no claim accounts for.
#      Runs after 01–10 in the same session. Output: output/tables/11_text_check.xlsx
# =============================================================================
source(file.path(DATA_DIR, "R_workflow_V25", "R", "read_docx_tables.R"))
P <- read_docx_paragraphs(MANUSCRIPT_DOCX)
i0 <- which(startsWith(P, "Summary"))[1]; i1 <- which(startsWith(P, "Contributors"))[1]
P <- P[i0:(i1 - 1)]
J <- jsonlite::read_json(file.path(OUT_DIR, "manuscript_values.json")); V <- J$V
d <- dat; bb <- filter(d, bbs == 1); nb <- filter(d, bbs == 0); ad <- filter(d, age18 %in% ">=18"); ch <- filter(d, age18 %in% "<18")
f1 <- function(x) md(x, 1); pc <- function(k, n) md(100 * k / n, 1); p0 <- function(k, n) md(100 * k / n, 0)
fish <- function(a, n1, b, n2) fisher.test(matrix(c(a, n1 - a, b, n2 - b), 2))$p.value
wci2 <- function(k, n) { w <- wilson(k, n); sprintf("%s–%s", md(w[2]), md(w[3])) }
s5 <- function(cmp, f, col = "p_adj_within") tableS5[[col]][tableS5$comparison == cmp & tableS5$feature == f]
C <- list(); add <- function(section, fragment, type = "R") C[[length(C) + 1]] <<- tibble(section = section, fragment = fragment, type = type)

# ---- Summary / Research in context -----------------------------------------------------------------
g <- sum(d$genetic); add("Summary", sprintf("%d participants (%s%%, 95%% CI %s) had a genetic disease", g, pc(g, 1004), wci2(g, 1004)))
add("Summary", sprintf("BBS in %d (%s%%), non-BBS obesity-associated disease in %d (%s%%), and non-obesity disease in %d (%s%%)",
    sum(d$bbs), pc(sum(d$bbs), 1004), sum(d$grp4 == 2), pc(sum(d$grp4 == 2), 1004), sum(d$grp4 == 3), pc(sum(d$grp4 == 3), 1004)))
nc <- table(factor(bb$bbs_cat, c("A", "B1", "B2", "C")))
add("Summary", sprintf("%d had biallelic pathogenic or likely pathogenic variants, %d one such variant", nc[["A"]], nc[["B1"]]))
add("Summary", sprintf("and %d a clinical assignment; %d came from two centres", nc[["B2"]] + nc[["C"]], sum(bb$high_site)))
add("Summary", sprintf("The Beales criteria identified %d of 38 (%s%%) and %d of 22 molecularly confirmed cases, the ERN criteria %d; %d had a prior clinical diagnosis",
    sum(bb$beales), pc(sum(bb$beales), 38), sum(bb$beales[bb$bbs_cat == "A"]), sum(bb$ern), sum(bb$clin_dx)))
add("Summary", sprintf("apparent area under the curve of %s (optimism-corrected %s)", V$auc_pt, V$auc_corr))
add("Summary", sprintf("at four or more points it identified %d of 38 cases (%s%%; specificity %s%%) and %d of 22", V$hi$tp, V$hi$sens_p, V$hi$spec_p, V$catA_hi))
add("Research in context", sprintf("%d (%s%%) had a genetic disease, %d (%s%%) an obesity-associated", g, pc(g, 1004), sum(d$obesity_gd), pc(sum(d$obesity_gd), 1004)))
add("Research in context", sprintf("Of 22 participants with molecularly confirmed BBS, %d did not meet the Beales criteria and %d did not meet the ERN",
    22 - sum(bb$beales[bb$bbs_cat == "A"]), 22 - sum(bb$ern[bb$bbs_cat == "A"])))
add("Research in context", sprintf("across all 38 study-assigned cases, %d had been diagnosed clinically", sum(bb$clin_dx)))
add("Research in context", sprintf("identified %d of 22 molecularly confirmed cases at a specificity of %s%%", V$catA_hi, V$hi$spec_p))

# ---- Methods ------------------------------------------------------------------------------------------
add("Methods", sprintf("not assessed in the %d participants", sum(!d$hyperphagia_assessed)))
add("Methods", sprintf("all comparisons of hyperphagia are restricted to the %d participants with the checklist", sum(d$hyperphagia_assessed)))
hb <- na.omit(bb$hyperphagia); hn <- na.omit(d$hyperphagia[d$grp4 == 3])
oth_ns <- all(as.numeric(gsub("<", "", chartr("·", ".", MT$table1[MT$table1[, 1] %in% c("Arterial hypertension", "BMI-SDS (AGA reference), children", "Age at enrolment, years, median (IQR)"), 7]))) >= 0.05)
add("Results", sprintf("Among participants with the checklist, hyperphagia was more frequent in BBS than in non-obesity genetic disease (%s%% vs %s%%; adjusted p=%s); arterial hypertension, BMI-SDS in children, and age did not differ between the three disease groups",
    pc(sum(hb), length(hb)), pc(sum(hn), length(hn)), md(tableS5$p_adj_within[tableS5$comparison == "BBS vs Non-obesity" & tableS5$feature == "hyperphagia"], 2)), if (oth_ns) "R" else "CHECK")
add("Methods", sprintf("maximum %d)", sum(SCORE_A_POINTS)))
add("Methods", sprintf("site-cluster bootstrap (2000 resamples of the %d sites)", n_distinct(d$site)))
add("Methods", sprintf("across %d features within each comparison", length(S5_VARS))); add("Methods", sprintf("across all %d tests", 6 * length(S5_VARS)))
add("Methods", "in 500 bootstrap resamples", "constant"); add("Methods", "(20 repeats)", "constant"); add("Methods", "with 16 candidate features", "constant")
add("Methods", "coefficients ≤0·25 not scored", "constant"); add("Methods", "validation at the other 35 sites", "constant")
add("Methods", "VUS-high (4–5 points), VUS-mid (2–3), and VUS-low (0–1)", "input"); add("Methods", "204 of 284 laboratory P/LP variants", "input")
add("Methods", "concordance for 30 laboratory P/LP primary findings", "input"); add("Methods", "upgrade of GALK1 p.(Pro28Thr) and three de-novo CNVs", "input")
add("Methods", "at 37 sites", "R"); add("Methods", "The BBS group comprised 38 participants", "R")

# ---- Results: participants ----------------------------------------------------------------------------
add("Results", "1009 individuals were enrolled; five were excluded (withdrawal of consent, n=1; missing biological material, n=1; sample receipt after study closure, n=3)", "input")
add("Results", sprintf("analytical cohort of %d participants (figure 1): %d index participants and %d relatives", nrow(d), sum(d$index), sum(d$index == 0)))
add("Results", "A further 289 relatives were tested", "input")
n18 <- sum(d$age18 %in% "<18"); n18b <- sum(d$age18 %in% ">=18"); n16 <- sum(d$age16 %in% "<16"); n16b <- sum(d$age16 %in% ">=16")
add("Results", sprintf("%d participants (%s%%) were younger than 18 years and %d (%s%%) 18 years or older; %d (%s%%) were younger than 16 years and %d (%s%%) 16 years or older",
    n18, pc(n18, 1004), n18b, pc(n18b, 1004), n16, pc(n16, 1004), n16b, pc(n16b, 1004)))
add("Results", sprintf("(age not recorded, n=%d; negative age in the eCRF, n=%d", sum(is.na(d$age_y)), sum(d$age_y < 0, na.rm = TRUE)))
q <- quantile(d$age_y, c(.5, .25, .75), na.rm = TRUE)
add("Results", sprintf("Median age was %s years (IQR %s–%s); %d (%s%%) participants were female and %d (%s%%)", f1(q[1]), f1(q[2]), f1(q[3]),
    sum(d$female), pc(sum(d$female), 1004), sum(d$white), pc(sum(d$white), 1004)))
add("Results", sprintf("Of %d adults, %d (%s%%) had a recorded ophthalmological feature (retinal dystrophy %d, visual loss only %d); for %d adults no ophthalmological feature was recorded, and a protocol deviation was documented for %d of them",
    nrow(ad), sum(ad$ophth), pc(sum(ad$ophth), nrow(ad)), sum(ad$retinal), sum(ad$visual & !ad$retinal), sum(!ad$ophth), sum(ad$pd[ad$ophth == 0])))
add("Results", sprintf("Of %d participants younger than 18 years, %d (%s%%) had a recorded ophthalmological feature and %d a protocol deviation",
    nrow(ch), sum(ch$ophth), pc(sum(ch$ophth), nrow(ch)), sum(ch$pd)))
ra <- filter(ad, route == "adult section only"); rb <- filter(ad, route == "both sections")
add("Results", sprintf("%d were documented with the adult eligibility section only", nrow(ra)))
add("Results", sprintf("%d with the paediatric checklist in addition (obesity before 6 years recorded in %d, protocol deviation in %d)", nrow(rb), sum(rb$ad_obese_lt6 %in% "yes"), sum(rb$pd)))
add("Results", sprintf("%d of the %d adults without a recorded ophthalmological feature belonged to this group", sum(rb$ophth == 0), sum(ad$ophth == 0)))

# ---- Results: frequency ---------------------------------------------------------------------------------
ci <- function(k) sprintf("%d (%s%%, %s)", k, pc(k, 1004), wci2(k, 1004))
add("Results", sprintf("%s BBS, %s a non-BBS obesity-associated genetic disease, and %s a non-obesity", ci(sum(d$bbs)), ci(sum(d$grp4 == 2)), ci(sum(d$grp4 == 3))))
add("Results", sprintf("%d (%s%%) had no qualifying primary genetic diagnosis", sum(d$grp4 == 4), pc(sum(d$grp4 == 4), 1004)))
add("Results", sprintf("present in %d participants (%s%%, %s)", sum(d$obesity_gd), pc(sum(d$obesity_gd), 1004), wci2(sum(d$obesity_gd), 1004)))
gk <- sum(ch$genetic); ga <- sum(ad$genetic)
add("Results", sprintf("genetic disease %d of %d, %s%%, vs %d of %d, %s%%, p=%s; BBS %s%% vs %s%%, p=%s", gk, nrow(ch), pc(gk, nrow(ch)), ga, nrow(ad), pc(ga, nrow(ad)),
    md(fish(gk, nrow(ch), ga, nrow(ad)), 2), pc(sum(ch$bbs), nrow(ch)), pc(sum(ad$bbs), nrow(ad)), md(fish(sum(ch$bbs), nrow(ch), sum(ad$bbs), nrow(ad)), 2)))
ao <- filter(ad, ophth == 1); an <- filter(ad, ophth == 0)
add("Results", sprintf("BBS was present in %d of %d (%s%%) with and %d of %d (%s%%) without a recorded ophthalmological feature (p=%s), and a genetic disease in %s%% and %s%%",
    sum(ao$bbs), nrow(ao), pc(sum(ao$bbs), nrow(ao)), sum(an$bbs), nrow(an), pc(sum(an$bbs), nrow(an)), md(fish(sum(ao$bbs), nrow(ao), sum(an$bbs), nrow(an)), 3),
    pc(sum(ao$genetic), nrow(ao)), pc(sum(an$genetic), nrow(an))))
add("Results", sprintf("BBS was present in %d of %d adults (%s%%) documented with the adult section only and in %d of %d (%s%%) documented with the paediatric checklist (p=%s)",
    sum(ra$bbs), nrow(ra), pc(sum(ra$bbs), nrow(ra)), sum(rb$bbs), nrow(rb), pc(sum(rb$bbs), nrow(rb)), md(fish(sum(ra$bbs), nrow(ra), sum(rb$bbs), nrow(rb)), 4)))
add("Results", if (all(rb$bbs_cat[rb$bbs == 1] == "C")) "both were in evidence category C" else "CHECK: category of the two paediatric-route adults")
add("Results", sprintf("%d of the %d came from one site", max(table(ra$site[ra$bbs == 1])), sum(ra$bbs)))
add("Results", "29 participants (2·9%) without a qualifying diagnosis carried an ACMG secondary finding (cancer predisposition 11, cardiomyopathy and arrhythmia 10, familial hypercholesterolaemia 7, transthyretin amyloidosis 1), similar to the 3·0%", "input")
hs <- filter(d, high_site == 1); os <- filter(d, high_site == 0)
add("Results", sprintf("Two of %d sites contributed %d of the 38 participants with BBS (%d participants enrolled at these sites); BBS was present in %d of %d participants (%s%%) at the other sites, and a genetic disease in %d (%s%%). BBS cases occurred at %d sites",
    n_distinct(d$site), sum(hs$bbs), nrow(hs), sum(os$bbs), nrow(os), pc(sum(os$bbs), nrow(os)), sum(os$genetic), pc(sum(os$genetic), nrow(os)), n_distinct(d$site[d$bbs == 1])))
add("Results", "contributed to seven diagnoses: five in BBS genes (BBS1, BBS7 in two participants, BBS9, BBS10) and two in group 3 (CFTR, EIF2B5); trans configuration was documented by segregation in three", "input")
add("Results", "Under a strict P/LP-only rule, 142 participants (14·1%)", "input"); add("Results", "128 (12·7%) had a genetic disease and 65 (6·5%)", "input")
add("Results", "BBS1 and BBS10 accounted for five each and MKKS for four; the recurrent alleles BBS10 c.271dup and BBS1 p.(Met390Arg) occurred in five and four participants", "input")
add("Results", "Among the 27 participants with BBS meeting the molecular diagnosis rule", "input")

add("Results", sprintf("Non-obesity genetic diseases (n=%d) comprised monogenic neurodevelopmental and neurological disorders in 24, CNVs in 12", sum(d$grp4 == 3)), "input")
add("Results", "metabolic, hearing, and other disorders in 27", "input"); add("Results", "CNVs or chromosomal aberrations caused 36 (25·0%) of the 144", "input")
add("Results", sprintf("Of the 38 participants with study-assigned BBS, %d had biallelic P/LP variants", nc[["A"]]))
add("Results", sprintf("%d a P/LP variant with a VUS-high second allele (B1), %d variants below the molecular rule with Beales or ERN criteria met (B2), and %d variants below the rule without either criterion (C)", nc[["B1"]], nc[["B2"]], nc[["C"]]))
add("Results", "phase not documented for three compound heterozygous genotypes", "input")
add("Results", "two VUS-mid or VUS-low alleles in eight participants, a homozygous VUS-high variant in two, and a P/LP variant with a VUS-mid variant in one", "input")

# ---- Results: criteria -----------------------------------------------------------------------------------
w <- wilson(sum(bb$beales), 38); bA <- filter(bb, bbs_cat == "A")
add("Results", sprintf("met by %d of 38 participants with study-assigned BBS (sensitivity %s%%, 95%% CI %s–%s) and by %d of 22 with molecularly confirmed BBS (%s%%); the ERN criteria by %d (%s%%) and %d (%s%%)",
    sum(bb$beales), pc(sum(bb$beales), 38), md(w[2]), md(w[3]), sum(bA$beales), pc(sum(bA$beales), 22), sum(bb$ern), pc(sum(bb$ern), 38), sum(bA$ern), pc(sum(bA$ern), 22)))
add("Results", sprintf("recorded before enrolment in %d participants (%s%%)", sum(bb$clin_dx), pc(sum(bb$clin_dx), 38)))
add("Results", sprintf("met by %d participants without BBS (positive predictive value %s%%): %d in group 2, %d in group 3, and %d in group 4",
    sum(nb$beales), pc(sum(bb$beales), sum(d$beales)), sum(d$beales[d$grp4 == 2]), sum(d$beales[d$grp4 == 3]), sum(d$beales[d$grp4 == 4])))
add("Results", if (sum(nb$ern) == 0) "The ERN criteria were met by no participant without BBS" else "CHECK ERN false positives")
add("Results", sprintf("the Beales criteria were met by %d and the ERN criteria by %d participants with BBS", sum(bb$beales_export), sum(bb$ern_export)))

# ---- Results: clinical characterisation ----------------------------------------------------------------------
qf <- quantile(bb$n_features, c(.5, .25, .75)); qd <- quantile(bb$n_domains, c(.5, .25, .75))
add("Results", sprintf("(median %d, IQR %d–%d) and affected domains (%d, %d–%d), with three or more domains in %d (%s%%)", qf[[1]], qf[[2]], qf[[3]], qd[[1]], qd[[2]], qd[[3]], sum(bb$dom3), pc(sum(bb$dom3), 38)))
add("Results", sprintf("Polydactyly (%s%%), retinal dystrophy (%s%%), and visual loss (%s%%)", pc(sum(bb$poly), 38), pc(sum(bb$retinal), 38), pc(sum(bb$visual), 38)))
B <- "BBS"; O <- "Non-BBS obesity-associated"; N <- "Non-obesity"; Q <- "No qualifying diagnosis"; cmpn <- function(a, b) paste(a, "vs", b)
sig <- function(a, b, f) s5(cmpn(a, b), f) < 0.05
add("Results", if (all(sapply(c("poly", "retinal", "visual"), function(f) sig(B, O, f) & sig(B, N, f) & sig(B, Q, f))))
  "were more frequent in BBS than in both other genetic groups and in group 4" else "CHECK poly/retinal/visual significance")
add("Results", sprintf("renal anomaly or dysfunction (%s%%) than in groups 2 and 4", pc(sum(bb$renal), 38)),
    if (sig(B, O, "renal") & sig(B, Q, "renal") & !sig(B, N, "renal")) "R" else "CHECK")
add("Results", sprintf("and hypogonadism (%s%%) and syndactyly (%s%%) than in groups 3 and 4", pc(sum(bb$hypog), 38), pc(sum(bb$syn), 38)),
    if (sig(B, N, "hypog") & sig(B, Q, "hypog") & sig(B, N, "syn") & sig(B, Q, "syn") & !sig(B, O, "hypog")) "R" else "CHECK")
add("Results", sprintf("(%s%% vs %s%%; adjusted p=%s)", pc(sum(bb$motor), 38), pc(sum(d$motor[d$grp4 == 3]), 63), md(s5(cmpn(B, N), "motor"), 2)))
g2 <- filter(d, grp4 == 2); g3_ <- filter(d, grp4 == 3); g4 <- filter(d, grp4 == 4)
orh <- function(a, n1, b, n2) md((a + .5) * (n2 - b + .5) / ((n1 - a + .5) * (b + .5)), 1)
add("Results", sprintf("predominantly male (%d of %d, %s%%)", sum(g2$male), nrow(g2), pc(sum(g2$male), nrow(g2))))
add("Results", sprintf("hypogonadism (%s%% vs %s%%; odds ratio %s), intellectual disability (%s%% vs %s%%; %s), cognitive impairment (%s%% vs %s%%)",
    pc(sum(g2$hypog), 43), pc(sum(g4$hypog), 860), orh(sum(g2$hypog), 43, sum(g4$hypog), 860), pc(sum(g2$id_dis), 43), pc(sum(g4$id_dis), 860),
    orh(sum(g2$id_dis), 43, sum(g4$id_dis), 860), pc(sum(g2$cog), 43), pc(sum(g4$cog), 860)))
mx <- max(sapply(c("hypog", "id_dis", "cog", "speech", "motor"), function(f) s5(cmpn(O, Q), f)))
add("Results", sprintf("(all adjusted p≤%s)", md(ceiling(mx * 1000) / 1000, 3)))
add("Results", if (sum(g2$renal) == 0) "renal anomalies were not recorded" else "CHECK renal group 2")
s2a <- filter(d, sub2 == "2a"); s2b <- filter(d, sub2 == "2b")
h2a <- na.omit(s2a$hyperphagia); h2b <- na.omit(s2b$hyperphagia)
add("Results", sprintf("Hyperphagia was recorded in %d of %d assessed participants in group 2a and %d of %d in group 2b (p=%s), and neurodevelopmental impairment in %d of %d and %d of %d (p=%s)",
    sum(h2a), length(h2a), sum(h2b), length(h2b), md(fish(sum(h2a), length(h2a), sum(h2b), length(h2b)), 2),
    sum(s2a$nd), nrow(s2a), sum(s2b$nd), nrow(s2b), md(fish(sum(s2a$nd), nrow(s2a), sum(s2b$nd), nrow(s2b)), 2)))
g2f <- g2_find %>% filter(id %in% s2a$id)
cnt_gene <- function(pat) n_distinct(g2f$id[str_detect(g2f$gene, pat)])
add("Results", sprintf("comprised %d disorders in which obesity is a core feature (group 2a: MC4R %d, LEPR %d, recurrent 16p11.2 deletions %d, 6q16 deletion including SIM1 %d, Smith–Magenis syndrome %d, PHIP %d) and %d syndromic disorders with variable obesity",
    nrow(s2a), cnt_gene("^MC4R$"), cnt_gene("^LEPR$"), cnt_gene("16p11\\.2\\("), cnt_gene("6q1[56]"), cnt_gene("17p11\\.2\\("), cnt_gene("^PHIP$"), nrow(s2b)))
add("Results", "(group 2b: PROK2 4, other CNVs and chromosomal disorders 11, monogenic syndromic or endocrine disorders 5)", "input")
add("Discussion", sprintf("accounted for %d of %d diagnoses", cnt_gene("^MC4R$") + cnt_gene("^LEPR$"), sum(d$grp4 == 2)))
add("Discussion", sprintf("heterozygous MC4R variants accounted for %s of %d", c("zero","one","two","three","four","five","six","seven","eight","nine","ten","eleven")[cnt_gene("^MC4R$") + 1], sum(d$grp4 == 2)))
mx3 <- max(sapply(c("cog", "id_dis", "speech", "motor", "ataxia"), function(f) s5(cmpn(N, Q), f)))
add("Results", sprintf("and ataxia (%s%% vs %s%%) than group 4 (adjusted p≤%s)", pc(sum(g3_$ataxia), 63), pc(sum(g4$ataxia), 860), md(ceiling(mx3 * 1000) / 1000, 3)))
add("Results", sprintf("dental anomalies were recorded in %s%% vs %s%% (adjusted p=%s)", pc(sum(g3_$dental), 63), pc(sum(g4$dental), 860), md(s5(cmpn(N, Q), "dental"), 2)))
add("Results", sprintf("(%s%% vs %s%%; adjusted p=%s within the comparison, %s across all tests)", pc(sum(g2$hypog), 43), pc(sum(g3_$hypog), 63),
    md(s5(cmpn(O, N), "hypog"), 2), md(s5(cmpn(O, N), "hypog", "p_adj_all144"), 3)))
add("Results", sprintf("all %d associations significant within comparisons remained significant", sum(tableS5$p_adj_within < .05)),
    if (all(tableS5$p_adj_all144[tableS5$p_adj_within < .05] < .05)) "R" else "CHECK")

# ---- Results: age strata ---------------------------------------------------------------------------------------
b1 <- filter(bb, age18 == "<18"); b2 <- filter(bb, age18 == ">=18")
add("Results", sprintf("retinal dystrophy was recorded in %d of %d (%s%%) participants younger than 18 years and %d of %d (%s%%) aged 18 years or older (p=%s), and visual loss in %s%% and %s%% (p=%s)",
    sum(b1$retinal), nrow(b1), pc(sum(b1$retinal), nrow(b1)), sum(b2$retinal), nrow(b2), pc(sum(b2$retinal), nrow(b2)),
    md(fish(sum(b1$retinal), nrow(b1), sum(b2$retinal), nrow(b2)), 2), p0(sum(b1$visual), nrow(b1)), p0(sum(b2$visual), nrow(b2)),
    md(fish(sum(b1$visual), nrow(b1), sum(b2$visual), nrow(b2)), 3)))
wt <- wilcox.test(age_y ~ retinal, data = bb, exact = FALSE)
add("Results", sprintf("(median %s vs %s years; p=%s)", f1(median(bb$age_y[bb$retinal == 0])), f1(median(bb$age_y[bb$retinal == 1])), md(wt$p.value, 3)))
add("Results", sprintf("Neurodevelopmental impairment (%s%% vs %s%%; p=%s), and intellectual disability (%s%% vs %s%%; p=%s)",
    p0(sum(b1$nd), nrow(b1)), p0(sum(b2$nd), nrow(b2)), md(fish(sum(b1$nd), nrow(b1), sum(b2$nd), nrow(b2)), 2),
    p0(sum(b1$id_dis), nrow(b1)), p0(sum(b2$id_dis), nrow(b2)), md(fish(sum(b1$id_dis), nrow(b1), sum(b2$id_dis), nrow(b2)), 2)))
add("Results", sprintf("assessed in %d of %d adults with BBS (present in %d)", sum(b2$hyperphagia_assessed), nrow(b2), sum(b2$hyperphagia, na.rm = TRUE)))
g4p <- sapply(c("retinal", "nd", "id_dis"), function(f) { a <- g4[g4$age18 %in% "<18", ][[f]]; b <- g4[g4$age18 %in% ">=18", ][[f]]
  fish(sum(a), length(a), sum(b), length(b)) })
add("Results", "Group 4 showed the same pattern (all p<0·0001)", if (all(g4p < 1e-4)) "R" else "CHECK")
add("Results", sprintf("%d of the %d adults with BBS and retinal dystrophy were documented with the adult eligibility section", sum(b2$retinal & b2$route == "adult section only"), sum(b2$retinal)))
na_ <- function(f) sum(d[[paste0(f, "_3")]] == "not_assessed")
add("Results", sprintf("retinal status was not assessed in %d participants, renal status in %d, genital status in %d, and polydactyly in %d", na_("retinal"), na_("renal"), na_("hypog_genital"), na_("poly")))
add("Results", sprintf("not assessed in %d of %d participants with BBS younger than 18 years and %d of %d adults", sum(b1$retinal_3 == "not_assessed"), nrow(b1), sum(b2$retinal_3 == "not_assessed"), nrow(b2)))

# ---- Results: score ---------------------------------------------------------------------------------------------
add("Results", sprintf("AUC of %s; the ordinal Beales feature count had an AUC of %s (difference %s, 95%% CI %s; DeLong p=%s; site-cluster 95%% CI %s",
    sub(" \\(", " (95% CI ", V$auc), V$auc_beales, V$diff, V$diff_ci, V$delong_p, V$cl_diff))
s14 <- setNames(J$S14, sapply(J$S14, `[[`, "item"))
add("Results", sprintf("polydactyly received a median of %s, retinal dystrophy %s, and hypogonadism %s; syndactyly or brachydactyly received at least one point in %s%%, neurodevelopmental impairment in %s%%, and dental anomaly in %s%%",
    sub(" \\(", " points (2·5th–97·5th percentile ", s14$poly$med), s14$retinal$med, s14$hypog_genital$med, s14$synbrachy$ge1, s14$nd$ge1, s14$dental$ge1))
add("Results", sprintf("Score A had an optimism-corrected AUC of %s; the full model had an optimism-corrected AUC of %s, a calibration slope of %s, a calibration-in-the-large of %s (%s after FLIC intercept correction), and a Brier score of %s (optimism-corrected %s)",
    V$auc_corr, V$auc_model_corr, V$slope, sub("^-", "−", md(validation$citl_firth_apparent, 2)), md(abs(validation$citl_flic_apparent), 2), V$brier, V$brier_corr))
add("Results", sprintf("corrected sensitivity of %s%% and specificity of %s%% at four or more points and of %s%% and %s%% at three or more points; a threshold re-chosen by Youden’s index in each resample was three to five points in %s%% of resamples (corrected sensitivity %s%%, specificity %s%%)",
    V$proc_sens_hi, V$proc_spec_hi, V$proc_sens_lo, V$proc_spec_lo, V$youden_pct_3to5, V$proc_sensY, V$proc_specY))
add("Results", sprintf("the AUC was %s (range %s), with a sensitivity of %s%% and a specificity of %s%% at four or more points", V$cv_auc, V$cv_range, V$cv_sens_hi, V$cv_spec_hi))
add("Results", sprintf("the score had an AUC of %s (95%% CI %s) and identified %s participants with BBS at four or more points (specificity %s%%)", V$site_auc, V$site_ci, sub("/", " of ", V$site_hi), V$site_spec))
add("Results", sprintf("Score A identified %d of 38 participants with BBS (sensitivity %s%%; specificity %s%%) and selected %d participants (%s%%); the observed probability of BBS below three points was %s%% (%s)",
    V$lo$tp, V$lo$sens_p, V$lo$spec_p, V$lo$sel, V$lo_sel_pct, V$below_lo_pct, V$below_lo))
add("Results", sprintf("it identified %d (sensitivity %s%%, 95%% CI %s; specificity %s%%, %s) and selected %d participants (%s%%), with a positive predictive value of %s%%",
    V$hi$tp, V$hi$sens_p, sub("^[^(]*\\((.*)\\)$", "\\1", V$hi$sens), V$hi$spec_p, sub("^[^(]*\\((.*)\\)$", "\\1", V$hi$spec), V$hi$sel, V$hi_sel_pct, V$hi$ppv_p))
add("Results", sprintf("predictive values would be %s%%, %s%%, and %s%%", V$ppv_prev[["2"]][[2]], V$ppv_prev[["1"]][[2]], V$ppv_prev[["0.5"]][[2]]))
add("Results", sprintf("the AUC was %s (95%% CI %s), the sensitivity at four or more points %s%% (%s), and the specificity %s%% (%s)", V$auc_pt, V$cl_auc, V$hi$sens_p, V$cl_sens_hi, V$hi$spec_p, V$cl_spec_hi))
add("Results", sprintf("met the Beales criteria and %d further participants (exact McNemar p=%s), with higher observed specificity (%s%% vs %s%%)",
    V$mcn_only_score, md(binom.test(V$mcn_only_score, V$mcn_only_score + V$mcn_only_beales)$p.value, 3), V$hi$spec_p, V$spec_beales))
add("Results", sprintf("Score A had an AUC of %s; at four or more points it identified %d of 22 (%s%%), compared with %d by the Beales and %d by the ERN",
    sub(" .*", "", V$catA_auc), V$catA_hi, pc(V$catA_hi, 22), sum(bA$beales), sum(bA$ern)))
add("Results", sprintf("the AUC was %s; leave-one-site-out AUCs ranged from %s", V$index_auc, V$loso))
add("Results", sprintf("Score A had an AUC of %s (95%% CI %s) and identified %d of 7 participants with BBS at four or more points, the Beales criteria %d",
    sub(" .*", "", V$oth_auc), sub("^[^(]*\\((.*)\\)$", "\\1", V$oth_auc), V$oth_hi, sum(os$beales[os$bbs == 1])))
add("Results", sprintf("identified %d of 23 participants with BBS (specificity %s%%), the Beales criteria %d, the ERN criteria %d, and a prior clinical diagnosis %d; in adults, %d, %d, %d, and %d of 15",
    V$lt18$tp, V$lt18$spec_p, sum(b1$beales), sum(b1$ern), sum(b1$clin_dx), V$ge18$tp, sum(b2$beales), sum(b2$ern), sum(b2$clin_dx)))
rs <- retinal_spec
add("Results", sprintf("specificity of %s%% in participants younger than 18 years and %s%% in adults", md(rs$spec_retinal[rs$age18 == "<18"]), md(rs$spec_retinal[rs$age18 == ">=18"])))
add("Results", sprintf("(%s%% in adults documented with the adult eligibility section, in whom it was an entry criterion, and %s%% in adults",
    md(100 * mean(ra$retinal[ra$bbs == 0] == 0), 0), md(100 * mean(rb$retinal[rb$bbs == 0] == 0))))
add("Results", sprintf("without it, Score A had an AUC of %s in adults", V$noret_adult_auc))
add("Results", sprintf("counted as present the AUC was %s, and in %d participants with all score items assessed %s", V$napres_auc, V$cc_n, V$cc_auc))
add("Results", sprintf("Among %d participants without BBS, a Firth model with 16 clinical features discriminated non-BBS genetic disease (groups 2–3; %d cases) with an optimism-corrected AUC of %s and non-BBS obesity-associated genetic disease (%d cases) with %s",
    nrow(nb), sum(nb$genetic), md(tableS6_model$auc_corrected[1], 2), sum(nb$grp4 == 2), md(tableS6_model$auc_corrected[2], 2)))
cnt <- nb$nd + nb$hypog_genital + nb$id_dis; s2_ <- cnt >= 2
add("Results", sprintf("of two or more selected %d participants (%s%%) and was associated with a probability of non-BBS genetic disease of %s%% (%s%% overall) and of non-BBS obesity-associated genetic disease of %s%% (%s%% overall), capturing %d of %d and %d of %d",
    sum(s2_), pc(sum(s2_), nrow(nb)), pc(sum(nb$genetic[s2_]), sum(s2_)), pc(sum(nb$genetic), nrow(nb)), pc(sum(nb$grp4[s2_] == 2), sum(s2_)),
    pc(sum(nb$grp4 == 2), nrow(nb)), sum(nb$genetic[s2_]), sum(nb$genetic), sum(nb$grp4[s2_] == 2), sum(nb$grp4 == 2)))

# ---- Discussion --------------------------------------------------------------------------------------------------
h4 <- na.omit(d$hyperphagia[d$grp4 == 4]); hbb <- na.omit(bb$hyperphagia); hnn <- na.omit(d$hyperphagia[d$grp4 == 3])
padj4 <- tableS5$p_adj_within[tableS5$comparison == "BBS vs No qualifying diagnosis" & tableS5$feature == "hyperphagia"]
add("Discussion", sprintf("hyperphagia was recorded in %s%% of those with BBS, more often than in non-obesity genetic disease (%s%%) but not significantly more often than in participants without a qualifying diagnosis (%s%%; adjusted p=%s)",
    pc(sum(hbb), length(hbb)), pc(sum(hnn), length(hnn)), pc(sum(h4), length(h4)), md(padj4, 2)), if (padj4 >= 0.05) "R" else "CHECK")
add("Discussion", sprintf("%s%% had a genetic disease: BBS in %s%%, a non-BBS obesity-associated genetic disease in %s%%, and a non-obesity genetic disease in %s%%",
    pc(g, 1004), pc(sum(d$bbs), 1004), pc(sum(d$grp4 == 2), 1004), pc(sum(d$grp4 == 3), 1004)))
add("Discussion", sprintf("did not identify %d of 22 participants with molecularly confirmed BBS and the ERN criteria, as operationalised, %d", 22 - sum(bA$beales), 22 - sum(bA$ern)))
add("Discussion", "four in five participants with BBS had not been diagnosed clinically", if (round((38 - sum(bb$clin_dx)) / 38, 1) == 0.8) "R" else "CHECK")
add("Discussion", "one in four was caused by a CNV", "input")
add("Discussion", sprintf("the frequency of BBS in adults (%s%%)", pc(sum(ad$bbs), nrow(ad))))
add("Discussion", sprintf("BBS was present in %s%% of participants with and %s%% without", pc(sum(ao$bbs), nrow(ao)), pc(sum(an$bbs), nrow(an))))
add("Discussion", sprintf("its specificity for BBS was %s%% in adults compared with %s%% in children; Score A without the retinal item had an AUC of %s in adults",
    md(rs$spec_retinal[rs$age18 == ">=18"]), md(rs$spec_retinal[rs$age18 == "<18"]), V$noret_adult_auc))
add("Discussion", sprintf("retinal dystrophy was not recorded in %d of %d participants with BBS", sum(b1$retinal == 0), nrow(b1)))
add("Discussion", sprintf("protocol deviation in %d of %d); in these adults BBS was present in %s%% and the specificity of the retinal item was %s%%, compared with %s%% and %s%%",
    sum(rb$pd), nrow(rb), pc(sum(rb$bbs), nrow(rb)), md(100 * mean(rb$retinal[rb$bbs == 0] == 0)), pc(sum(ra$bbs), nrow(ra)), md(100 * mean(ra$retinal[ra$bbs == 0] == 0), 0)))
add("Discussion", "no participant with BBS was aged 16 or 17 years", if (!any(bb$age_y >= 16 & bb$age_y < 18)) "R" else "CHECK")
add("Discussion", sprintf("In WEGIO, %d of 22 molecularly confirmed cases did not meet the Beales criteria", 22 - sum(bA$beales)))
add("Discussion", sprintf("(%d and %d participants with BBS when derived from the eligibility checklist alone)", sum(bb$beales_export), sum(bb$ern_export)))
add("Discussion", sprintf("The Beales criteria were in addition met by %d participants without BBS", sum(nb$beales)))
add("Discussion", sprintf("%d of the %d participants with polydactyly had BBS", V$n_poly_bbs, V$n_poly))
add("Discussion", sprintf("it identified %d of 38 participants with study-assigned BBS and %d of 22 with molecularly confirmed BBS, every participant identified by the Beales criteria and %d additional participants",
    V$hi$tp, V$catA_hi, V$mcn_only_score))
add("Discussion", sprintf("Below three points, the observed probability of BBS was %s%%", V$below_lo_pct))
add("Discussion", sprintf("%s participants with BBS scored below three points", c("zero", "one", "two", "three", "four", "five")[sum(bb$scoreA < 3) + 1]))
add("Discussion", sprintf("applied to the others (AUC %s)", V$site_auc))
add("Discussion", sprintf("The positive predictive value of %s%% applies to a population with a BBS prevalence of about 4%%; at 1%% it would be about %s%%", V$hi$ppv_p, V$ppv_prev[["1"]][[2]]))
add("Discussion", sprintf("BBS was present in %s%% of participants and the discrimination of Score A was lower and imprecise (AUC %s; seven cases)", pc(sum(os$bbs), nrow(os)), sub(" .*", "", V$oth_auc)))
add("Discussion", sprintf("the 95%% CI of the overall AUC extended from %s", sub("–", " to ", V$cl_auc)))
add("Discussion", sprintf("%d of 38 participants with study-assigned BBS had variants below the molecular rule", nc[["B2"]] + nc[["C"]]))
add("Discussion", "eight of them with two VUS-mid or VUS-low alleles and two with a homozygous VUS-high variant", "input")
add("Discussion", sprintf("%s of these met the Beales or ERN criteria, and %s met neither", c("zero", "one", "two", "three", "four", "five", "six")[nc[["B2"]] + 1], c("zero", "one", "two", "three", "four", "five")[nc[["C"]] + 1]))

add("Discussion", "ADCY3 (16 participants), POMC (ten), and MRAP2 (two)", "input")
add("Discussion", sprintf("(optimism-corrected AUC %s)", md(tableS6_model$auc_corrected[1], 2)))
add("Discussion", sprintf("%d relatives were included", sum(d$index == 0)))
add("Discussion", sprintf("With %d BBS cases, of which %d were molecularly confirmed", sum(d$bbs), nc[["A"]]))


# ---- additional data claims (numbers repeated in several places) ------------------------------------
add("Summary", sprintf("In this multicentre cohort study, %d participants", nrow(d))); add("Research in context", sprintf("case-report forms in %d participants", nrow(d)))
add("Discussion", sprintf("In this multicentre cohort of %d participants", nrow(d)))
add("Methods", sprintf("the sensitivity of the criteria against all %d participants", sum(d$bbs)))
add("Methods", sprintf("assigns polydactyly %d points, retinal dystrophy %d, and hypogonadism or genital hypoplasia, renal anomaly or dysfunction, syndactyly or brachydactyly, and neurodevelopmental impairment %d each (maximum %d)",
    SCORE_A_POINTS[["poly"]], SCORE_A_POINTS[["retinal"]], SCORE_A_POINTS[["renal"]], sum(SCORE_A_POINTS)))
add("Methods", sprintf("performed in R version %s.%s", R.version$major, R.version$minor))
add("Results", sprintf("The Firth model gave %d points for polydactyly, %d for retinal dystrophy, and %d each", SCORE_A_POINTS[["poly"]], SCORE_A_POINTS[["retinal"]], SCORE_A_POINTS[["nd"]]))
add("Results", sprintf("received %d instead of %d points. In %d bootstrap resamples", SCORE_A_POINTS[["renal"]], SCORE_A_LEGACY[["renal"]], N_BOOT))
add("Results", sprintf("applied to the other %d sites (%d BBS)", n_distinct(d$site) - 2, sum(os$bbs)))
add("Results", sprintf("Because %d of the %d participants with BBS came from two sites", sum(hs$bbs), sum(d$bbs)))
add("Results", sprintf("At the %d sites outside the two sites", n_distinct(d$site) - 2))
add("Results", sprintf("With correction across all %d pairwise tests", 6 * length(S5_VARS)))
add("Results", sprintf("Non-BBS obesity-associated genetic diseases (n=%d) comprised", sum(d$grp4 == 2)))
add("Results", sprintf("Of the %d study-assigned BBS cases", sum(d$bbs)))
add("Results", "at prevalences of 2%, 1%, and 0·5%", "constant"); add("Results", "at the cohort prevalence of 3·8%", "R")
add("Discussion", sprintf("Two sites contributed %d of the %d participants with BBS. At the other %d sites", sum(hs$bbs), sum(d$bbs), n_distinct(d$site) - 2))

claims <- bind_rows(C)
body <- paste(P, collapse = "\n")
claims <- claims %>% mutate(found = vapply(fragment, function(f) grepl(f, body, fixed = TRUE), logical(1)),
  status = case_when(grepl("^CHECK", fragment) ~ "R CHECK FAILED", type == "R" & found ~ "match", type == "R" ~ "MISMATCH",
                     type == "CHECK" ~ "R CHECK FAILED", found ~ paste(type, "(present)"), TRUE ~ paste(type, "(NOT FOUND)")))

# ---- coverage: numbers in the checked sections that no claim accounts for --------------------------------------
sec <- P
keep <- !grepl("^(Introduction|Obesity affects|Genetic obesity|Genetic testing nevertheless|BBS illustrates|We searched PubMed|Funding|WEGIO is a multicentre|A family-member|Population data support)", sec)
num_re <- "(?<![A-Za-z0-9·])[0-9]+(?:·[0-9]+)?(?![A-Za-z0-9])"
cover <- map_dfr(which(keep), function(i) { p <- sec[i]
  hits <- gregexpr(num_re, p, perl = TRUE)[[1]]; if (hits[1] < 0) return(NULL)
  pos <- as.integer(hits); len <- attr(hits, "match.length")
  covered <- logical(nchar(p))
  for (f in claims$fragment[claims$found]) { m <- gregexpr(f, p, fixed = TRUE)[[1]]
    if (m[1] > 0) for (k in seq_along(m)) covered[m[k]:(m[k] + nchar(f) - 1)] <- TRUE }
  tibble(paragraph = i, number = substring(p, pos, pos + len - 1), covered = covered[pos],
         context = substring(p, pmax(1, pos - 60), pmin(nchar(p), pos + len + 40)),
         pre = substring(p, pmax(1, pos - 40), pos - 1), post = substring(p, pos + len, pmin(nchar(p), pos + len + 20)))
}) %>% filter(!covered)
classify <- function(pre, post, num) {
  pre <- sub("\\s+$", "", pre, perl = TRUE); pre_raw <- pre
  if (grepl("(?i)(table|tables|figure|appendix|group|groups|category|S)\\s*$", pre_raw, perl = TRUE) ||
      grepl("(?i)(table|figure|appendix) S?[0-9]+(,| and) (table )?S?$", pre_raw, perl = TRUE)) return("table/figure/appendix/group reference")
  if (grepl("[A-Za-z.;,)]$", pre) && !grepl("\\s$", pre)) return("reference citation")
  if (grepl("^[–,][0-9]", post) && grepl("[A-Za-z.;,)][0-9]*$", pre)) return("reference citation")
  if (grepl("\\($", pre) && grepl("^\\)", post)) return("enumeration")
  if (grepl("^[a-z]\\b|^[a-z]:", post)) return("subgroup label (e.g. 2a)")
  if (grepl("(groups? [0-9]+ and|groups [0-9]+–|[0-9]–)$", pre) ) return("table/figure/appendix/group reference")
  if (grepl("V[0-9]*·?$|V$", pre) || grepl("^·[0-9]", post) && grepl("V$", pre)) return("version number")
  if (grepl("^( years| years|th percentile| kg/m|-year| ye)", post) || grepl("(<|≥|BMI of|splits at|split at|at|than|from|before)$", sub("\\s+$", "", pre)))
    return("definition (eligibility / age split)")
  if (num == "95") return("confidence level")
  "DATA NUMBER NOT COVERED"
}
cover$reason <- if (nrow(cover)) mapply(classify, cover$pre, cover$post, cover$number) else character(0)

write_xlsx(list(claims = claims, uncovered_numbers = cover), file.path(OUT_DIR, "tables", "11_text_check.xlsx"))
message("11 text check: ", nrow(claims), " claims; ", sum(claims$status == "match"), " match; ",
        sum(claims$status %in% c("MISMATCH", "R CHECK FAILED")), " mismatch/failed; ",
        sum(grepl("NOT FOUND", claims$status)), " input/constant not found; ", sum(cover$reason == "DATA NUMBER NOT COVERED"),
        " data numbers not covered (", nrow(cover), " structural numbers classified)")
