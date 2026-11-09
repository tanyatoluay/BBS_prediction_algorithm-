# =============================================================================
# 03 — Clinical criteria and binary rules: Table 4, S2, S11, McNemar
# =============================================================================
dat <- readRDS(file.path(OUT_DIR, "wegio_analysis_dataset.rds"))

# diagnostic accuracy of a binary rule `test` against reference `ref`
diag_acc <- function(test, ref) {
  tp <- sum(test == 1 & ref == 1); fn <- sum(test == 0 & ref == 1)
  fp <- sum(test == 1 & ref == 0); tn <- sum(test == 0 & ref == 0)
  se <- wilson(tp, tp + fn); sp <- wilson(tn, tn + fp); ppv <- wilson(tp, tp + fp)
  tibble(identified = sprintf("%d/%d", tp, tp + fn),
         sens = sprintf("%.1f (%.1f–%.1f)", se[1], se[2], se[3]),
         spec = sprintf("%.1f (%.1f–%.1f)", sp[1], sp[2], sp[3]),
         selected = tp + fp, ppv = round(ppv[1], 1), npv = round(100 * tn / (tn + fn), 1),
         tp = tp, fn = fn, fp = fp, tn = tn)
}
RULES  <- c(clin_dx = "Clinical diagnosis before enrolment", ern = "ERN criteria", beales = "Beales criteria",
            beales_v25 = "Beales (V25 counting)", scoreA3 = "Score A ≥3", scoreA5 = "Score A ≥5")
STRATA <- list(All = rep(TRUE, nrow(dat)), `<18` = dat$age18 %in% "<18", `>=18` = dat$age18 %in% ">=18",
               `<16` = dat$age16 %in% "<16", `>=16` = dat$age16 %in% ">=16")
table4 <- imap_dfr(RULES, function(lab, r) imap_dfr(STRATA, function(s, sl)
  diag_acc(dat[[r]][s], dat$bbs[s]) %>% mutate(criterion = lab, stratum = sl, .before = 1)))

# McNemar (exact) among BBS: Score A>=5 vs Beales
mcn <- function(a, b) { tb <- table(factor(a, 0:1), factor(b, 0:1));
  list(b10 = tb[2, 1], b01 = tb[1, 2], p = binom.test(tb[2, 1], tb[2, 1] + tb[1, 2])$p.value) }
bb <- filter(dat, bbs == 1)
mcnemar_res <- tibble(
  comparison = c("Score A>=5 vs Beales", "Score A>=5 vs ERN", "Beales vs ERN"),
  only_first = c(mcn(bb$scoreA5, bb$beales)$b10, mcn(bb$scoreA5, bb$ern)$b10, mcn(bb$beales, bb$ern)$b10),
  only_second = c(mcn(bb$scoreA5, bb$beales)$b01, mcn(bb$scoreA5, bb$ern)$b01, mcn(bb$beales, bb$ern)$b01),
  exact_p = c(mcn(bb$scoreA5, bb$beales)$p, mcn(bb$scoreA5, bb$ern)$p, mcn(bb$beales, bb$ern)$p))

# Table S11: evidence categories
tableS11 <- bb %>% group_by(bbs_cat) %>% summarise(n = n(), lt18 = sum(age18 == "<18"), beales = sum(beales),
  ern = sum(ern), scoreA5 = sum(scoreA5), clin_dx = sum(clin_dx), ids = paste(id, collapse = ", "))

# Table S2: export (eligibility checklist only) vs all items
s2 <- function(v) c(all = sum(dat[[v]]), bbs = sum(dat[[v]][dat$bbs == 1]),
  bbs_lt16 = sum(dat[[v]][dat$bbs == 1 & dat$age16 %in% "<16"]), bbs_ge16 = sum(dat[[v]][dat$bbs == 1 & dat$age16 %in% ">=16"]),
  nonbbs = sum(dat[[v]][dat$bbs == 0]), grp4 = sum(dat[[v]][dat$grp4 == 4]),
  spec = round(100 * mean(dat[[v]][dat$bbs == 0] == 0), 1))
tableS2 <- bind_rows(beales_export = s2("beales_export"), beales_all = s2("beales"), beales_v25 = s2("beales_v25"),
                     ern_export = s2("ern_export"), ern_all = s2("ern"), .id = "rule")
export_not_all <- dat %>% filter(beales_export == 1, beales == 0) %>% select(id, group, beales_prim, beales_sec, exp_beales_prim_sym)

# criteria by molecular status (table 4 footnote)
by_molecular <- bb %>% mutate(mol = if_else(bbs_cat %in% c("A", "B1"), "rule met", "below rule")) %>%
  group_by(mol) %>% summarise(n = n(), beales = sum(beales), scoreA3 = sum(scoreA3), scoreA5 = sum(scoreA5))

write_xlsx(list(table4 = table4, mcnemar = mcnemar_res, tableS11 = tableS11, tableS2 = tableS2,
                export_not_all_items = export_not_all, by_molecular = by_molecular),
           file.path(OUT_DIR, "tables", "03_criteria.xlsx"))
message("03 criteria done")
