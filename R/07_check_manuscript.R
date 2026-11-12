# =============================================================================
# 07 (checks the LEGACY Score A used in V25; the re-derived score is the new primary)
# 07 — Automated consistency check: V25 manuscript values vs values reproduced here
#      (run after 01–05 in the same session; see run_all.R)
# =============================================================================
d <- readRDS(file.path(OUT_DIR, "wegio_analysis_dataset.rds"))
ad <- filter(d, age18 == ">=18"); bb <- filter(d, bbs == 1); nb <- filter(d, bbs == 0)
r1 <- function(x) round(x, 1); r3 <- function(x) round(x, 3)
sens <- function(v, s = TRUE) sum(d[[v]][d$bbs == 1 & s])
spec <- function(v, s = TRUE) r1(100 * mean(d[[v]][d$bbs == 0 & s] == 0))
gA <- d %>% filter(is.na(bbs_cat) | bbs_cat == "A") %>% mutate(bbs = as.integer(bbs_cat %in% "A"))

chk <- tribble(
  ~section, ~item, ~manuscript, ~reproduced,
  "Cohort", "N analysed", 1004, nrow(d),
  "Cohort", "index participants", 983, sum(d$index),
  "Cohort", "<18 years", 670, sum(d$age18 == "<18", na.rm = TRUE),
  "Cohort", ">=18 years", 333, sum(d$age18 == ">=18", na.rm = TRUE),
  "Cohort", "<16 years", 623, sum(d$age16 == "<16", na.rm = TRUE),
  "Cohort", "female", 519, sum(d$female),
  "Cohort", "White", 877, sum(d$white),
  "Cohort", "median age", 12.6, r1(median(d$age_y, na.rm = TRUE)),
  "Yield", "BBS", 38, sum(d$grp4 == 1),
  "Yield", "non-BBS obesity", 43, sum(d$grp4 == 2),
  "Yield", "non-obesity", 63, sum(d$grp4 == 3),
  "Yield", "genetic disease", 144, sum(d$genetic),
  "Yield", "genetic disease <18", 98, sum(d$genetic[d$age18 %in% "<18"]),
  "Yield", "genetic disease >=18", 46, sum(d$genetic[d$age18 %in% ">=18"]),
  "Ophth", "adults with ophth feature", 165, sum(ad$ophth),
  "Ophth", "adults without ophth feature", 168, sum(!ad$ophth),
  "Ophth", "  of whom protocol deviation", 166, sum(ad$pd[ad$ophth == 0]),
  "Ophth", "children with ophth feature", 60, sum(d$ophth[d$age18 %in% "<18"]),
  "Ophth", "children with protocol deviation", 19, sum(d$pd[d$age18 %in% "<18"]),
  "Ophth", "BBS % adults with ophth", 7.9, r1(100 * mean(ad$bbs[ad$ophth == 1])),
  "Ophth", "BBS % adults without ophth", 1.2, r1(100 * mean(ad$bbs[ad$ophth == 0])),
  "Ophth", "p (Fisher)", 0.003, r3(fisher.test(table(ad$ophth, ad$bbs))$p.value),
  "Ophth", "retinal specificity <18", 96.8, spec("retinal", d$age18 %in% "<18"),
  "Ophth", "retinal specificity >=18", 58.2, spec("retinal", d$age18 %in% ">=18"),
  "Criteria", "Beales met, all", 89, sum(d$beales),
  "Criteria", "Beales met, BBS", 26, sens("beales"),
  "Criteria", "Beales met, category A", 16, sum(bb$beales[bb$bbs_cat == "A"]),
  "Criteria", "Beales met, non-BBS", 63, sum(nb$beales),
  "Criteria", "Beales met, group 4", 47, sum(d$beales[d$grp4 == 4]),
  "Criteria", "Beales PPV %", 29.2, r1(100 * mean(d$bbs[d$beales == 1])),
  "Criteria", "Beales specificity %", 93.5, spec("beales"),
  "Criteria", "Beales (export) met, BBS", 19, sens("beales_export"),
  "Criteria", "Beales (export) met, all", 35, sum(d$beales_export),
  "Criteria", "ERN met, BBS", 15, sens("ern"),
  "Criteria", "ERN met, category A", 11, sum(bb$ern[bb$bbs_cat == "A"]),
  "Criteria", "ERN (export) met, BBS", 4, sens("ern_export"),
  "Criteria", "prior clinical diagnosis, BBS", 7, sens("clin_dx"),
  "Score A", "AUC all BBS", 0.941, r3(auc_fast(d$scoreA_legacy, d$bbs)),
  "Score A", "AUC category A (others excluded)", 0.996, r3(auc_fast(gA$scoreA_legacy, gA$bbs)),
  "Score A", "AUC Beales count", 0.908, r3(auc_fast(d$beales_count, d$bbs)),
  "Score A", ">=5 identified", 33, sens("scoreA_legacy5"),
  "Score A", ">=5 specificity", 97.5, spec("scoreA_legacy5"),
  "Score A", ">=5 category A", 21, sum(bb$scoreA_legacy5[bb$bbs_cat == "A"]),
  "Score A", ">=3 identified", 34, sens("scoreA_legacy3"),
  "Score A", ">=3 specificity", 83.2, spec("scoreA_legacy3"),
  "Score A", "BBS below 3 points (of 808)", 4, sum(d$bbs[d$scoreA_legacy < 3]),
  "Score A", "AUC index only", 0.940, r3(with(filter(d, index == 1), auc_fast(scoreA_legacy, bbs))),
  "Score A", "AUC outside sites 10/39", 0.769, r3(with(filter(d, high_site == 0), auc_fast(scoreA_legacy, bbs))),
  "Score A", "BBS outside sites 10/39", 7, sum(d$bbs[d$high_site == 0]),
  "Score A", "N outside sites 10/39", 935, sum(d$high_site == 0),
  "Score A", "AUC >=18 retinal removed", 0.945, r3(with(filter(d, age18 %in% ">=18"), auc_fast(scoreA_legacy - 2 * retinal, bbs))),
  "Score A", "AUC not-assessed as present", 0.916, r3(auc_fast(with(d, 5*(poly_3!="absent")+2*(retinal_3!="absent")+2*(hypog_genital_3!="absent")+2*(renal_3!="absent")+(synbrachy_3!="absent")+(nd_3!="absent")), d$bbs)),
  "Score A", "complete-case N", 619, sum(d$scoreA_complete),
  "Score A", "complete-case AUC", 0.945, r3(with(filter(d, scoreA_complete), auc_fast(scoreA_legacy, bbs))),
  "Model", "Firth beta polydactyly", 5.14, round(unname(coef(fit)["poly"]), 2),
  "Model", "Firth beta renal", 1.38, round(unname(coef(fit)["renal"]), 2),
  "Model", "optimism-corrected AUC, full model", 0.933, r3(validation$auc_model_corrected),
  "Model", "optimism-corrected AUC, re-derived score", 0.930, r3(validation$auc_score_rederived_corrected),
  "Model", "calibration slope (corrected)", 0.97, round(validation$cal_slope_corrected, 2),
  "Model", "Brier", 0.012, r3(validation$brier_apparent),
  "Model", "bootstrap median points polydactyly", 6, tableS14$median[tableS14$item == "poly"],
  "Model", "bootstrap median points hypogonadism", 2, tableS14$median[tableS14$item == "hypog_genital"],
  "Model", "DeLong p", 0.0001, signif(roc.test(roc(d$bbs, d$scoreA_legacy, direction = "<", quiet = TRUE), roc(d$bbs, d$beales_count, direction = "<", quiet = TRUE), method = "delong", paired = TRUE)$p.value, 2),
  "Model", "LOSO AUC min", 0.911, r3(min(sapply(unique(d$site), function(s) with(d[d$site != s, ], auc_fast(scoreA_legacy, bbs))))),
  "Model", "LOSO AUC max", 0.961, r3(max(sapply(unique(d$site), function(s) with(d[d$site != s, ], auc_fast(scoreA_legacy, bbs))))),
  "Non-BBS", "count>=2 selected", 163, tableS6_count$n[tableS6_count$count_ge == 2],
  "Non-BBS", "optimism-corrected AUC non-BBS GD", 0.64, round(tableS6_model$auc_corrected[1], 2)
) %>% mutate(match = case_when(abs(manuscript - reproduced) < 1e-9 ~ "exact",
                               abs(manuscript - reproduced) <= 0.0051 * pmax(1, abs(manuscript)) ~ "rounding/MC",
                               TRUE ~ "DIFFERS"))
write_xlsx(chk, file.path(OUT_DIR, "tables", "07_manuscript_consistency_check.xlsx"))
if (!isTRUE(getOption("knitr.in.progress"))) print(as.data.frame(chk), right = FALSE)
message("07 check: ", sum(chk$match == "DIFFERS"), " of ", nrow(chk), " items differ")
