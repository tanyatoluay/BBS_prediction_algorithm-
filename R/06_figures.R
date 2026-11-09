# =============================================================================
# 06 — Figures: figure 2, figure 3 (A ROC, B by eligibility stratum), figure S1 calibration,
#      new figure S2 (adults by documented eligibility route)
# =============================================================================
dat  <- readRDS(file.path(OUT_DIR, "wegio_analysis_dataset.rds"))
obj  <- readRDS(file.path(OUT_DIR, "scoreA_objects.rds"))
cal  <- readRDS(file.path(OUT_DIR, "calibration.rds"))
FIG  <- file.path(OUT_DIR, "figures")
theme_w <- theme_classic(base_size = 10) + theme(legend.position = "top", legend.title = element_blank())
grp_fill <- c("BBS" = "black", "Non-BBS obesity-associated" = "grey40", "Non-obesity" = "white",
              "No qualifying diagnosis" = "grey82")

# ---- Figure 2: clinical features by genetic group --------------------------------------------------------
F2 <- c(poly = "Polydactyly", retinal = "Retinal dystrophy", hypog_genital = "Hypogonadism or micropenis",
        renal = "Renal anomaly or dysfunction", synbrachy = "Syndactyly or brachydactyly",
        motor = "Delayed motor development", cog = "Cognitive impairment", id_dis = "Intellectual disability",
        speech = "Speech delay or disorder", hyperphagia = "Hyperphagia", dental = "Dental anomaly",
        ataxia = "Ataxia", dom3 = "≥3 organ-system domains")
f2 <- dat %>% select(group, all_of(names(F2))) %>% pivot_longer(-group) %>%
  group_by(group, name) %>% summarise(pct = 100 * mean(value), n = n(), .groups = "drop") %>%
  mutate(feature = factor(F2[name], levels = rev(F2)),
         group_lab = factor(sprintf("%s (n=%d)", group, n), levels = sprintf("%s (n=%d)", levels(dat$group), table(dat$group))))
p2 <- ggplot(f2, aes(pct, feature, fill = group_lab)) +
  geom_col(position = position_dodge2(reverse = TRUE, width = .85), width = .85, colour = "black", linewidth = .2) +
  scale_fill_manual(values = setNames(unname(grp_fill), levels(f2$group_lab))) +
  scale_x_continuous(limits = c(0, 100), expand = c(0, 0)) + guides(fill = guide_legend(nrow = 2)) +
  labs(x = "Participants with the feature recorded (%)", y = NULL) + theme_w
ggsave(file.path(FIG, "Figure2_clinical_profile.png"), p2, width = 8, height = 7.5, dpi = 300)

# ---- Figure 3A: ROC with operating points ---------------------------------------------------------------------
roc_df <- function(r, lab) tibble(fpr = 1 - r$specificities, tpr = r$sensitivities, curve = lab) %>% arrange(fpr, tpr)
rocs <- bind_rows(roc_df(obj$rocA, sprintf("Score A (AUC %.3f)", auc(obj$rocA))),
                  roc_df(obj$rocB, sprintf("Beales feature count (AUC %.3f)", auc(obj$rocB))))
op <- function(v, lab) tibble(fpr = mean(dat[[v]][dat$bbs == 0]), tpr = mean(dat[[v]][dat$bbs == 1]), rule = lab)
ops <- bind_rows(op("beales", "Beales criteria met"), op("ern", "ERN criteria met"),
                 op("clin_dx", "Clinical diagnosis before enrolment"))
thr <- bind_rows(op("scoreA5", "≥5"), op("scoreA3", "≥3"))
p3a <- ggplot(rocs, aes(fpr, tpr)) + geom_abline(colour = "grey75", linewidth = .3) +
  geom_step(aes(linetype = curve), direction = "vh", linewidth = .7) +
  geom_point(data = thr, size = 2.4) + geom_text(data = thr, aes(label = rule), nudge_x = .05, nudge_y = -.04, size = 3) +
  geom_point(data = ops, aes(shape = rule), size = 3, fill = "white") +
  scale_shape_manual(values = c(23, 24, 22)) + coord_equal() +
  labs(x = "1 − specificity", y = "Sensitivity",
       title = sprintf("A  Discrimination of BBS (N=%d; %d BBS)", nrow(dat), sum(dat$bbs))) +
  theme_w + theme(legend.position = c(.62, .25), legend.direction = "vertical", legend.background = element_blank())

# ---- Figure 3B: sensitivity by eligibility stratum ---------------------------------------------------------------------
RULES_B <- c(clin_dx = "Prior clinical\ndiagnosis", ern = "ERN\ncriteria", beales = "Beales\ncriteria",
             scoreA5 = "Score A\n≥5", scoreA3 = "Score A\n≥3")
f3b <- dat %>% filter(bbs == 1) %>% select(age18, all_of(names(RULES_B))) %>% pivot_longer(-age18) %>%
  group_by(age18, name) %>% summarise(k = sum(value), n = n(), .groups = "drop") %>%
  mutate(pct = 100 * k / n, rule = factor(RULES_B[name], levels = RULES_B),
         stratum = sprintf("%s years (n=%d)", ifelse(age18 == "<18", "<18", "≥18"), n))
p3b <- ggplot(f3b, aes(rule, pct, fill = stratum)) +
  geom_col(position = position_dodge(.8), width = .8, colour = "black", linewidth = .2) +
  geom_text(aes(label = k), position = position_dodge(.8), vjust = -.4, size = 3) +
  scale_fill_manual(values = c("black", "grey75")) + scale_y_continuous(limits = c(0, 105), expand = c(0, 0)) +
  labs(x = NULL, y = "BBS participants identified (%)", title = "B  Sensitivity by eligibility stratum") +
  theme_w + theme(legend.position = c(.3, .92))
ggsave(file.path(FIG, "Figure3_scoreA_criteria.png"), p3a + p3b, width = 11, height = 5.2, dpi = 300)

# ---- Figure S1: calibration — (a) V25 decile plot, (b) recommended smoothed curve, Firth vs FLIC ----------------------
s1a <- ggplot(cal$calib, aes(pred, obs)) + geom_abline(colour = "grey60") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = .006) + geom_point(size = 2) + coord_equal(xlim = c(0, .6), ylim = c(0, .6)) +
  labs(x = "Predicted probability of BBS (decile mean)", y = "Observed proportion with BBS", title = "a  Deciles (as in V25)") + theme_w
sm <- bind_rows(tibble(p = cal$p, y = cal$y, model = "Firth"), tibble(p = cal$p_flic, y = cal$y, model = "Firth + FLIC intercept"))
s1b <- ggplot(sm, aes(pmax(p, 1e-4), y, linetype = model)) +
  geom_line(data = tibble(x = 10^seq(-4, 0, length.out = 200)), aes(x, x), inherit.aes = FALSE, colour = "grey60") +
  geom_smooth(method = "loess", se = FALSE, span = .9, colour = "black", formula = y ~ x) +
  geom_rug(data = filter(sm, model == "Firth", y == 0), sides = "b", alpha = .2) +
  geom_rug(data = filter(sm, model == "Firth", y == 1), sides = "t", alpha = .6) +
  scale_x_log10("Predicted probability (log scale)", breaks = c(.001, .01, .05, .2, .5, 1), labels = c(".001", ".01", ".05", ".2", ".5", "1")) +
  scale_y_continuous("Observed proportion (loess)") + coord_cartesian(ylim = c(0, 1)) +
  labs(title = "b  Smoothed calibration (recommended)") + theme_w
# loess on probability scale vs predicted: draw on probability y, logit x
ggsave(file.path(FIG, "FigureS1_calibration.png"), s1a + s1b, width = 11, height = 5, dpi = 300)

# ---- New figure S2: adults by documented eligibility route ----------------------------------------------------------
rt <- dat %>% filter(age18 == ">=18", route %in% c("adult section only", "both sections")) %>%
  mutate(route = recode(route, "adult section only" = "Adult criteria\n(adult section only)",
                        "both sections" = "Paediatric criteria\n(both sections; protocol deviation)")) %>%
  group_by(route) %>% summarise(n = n(), `Retinal dystrophy` = 100 * mean(retinal), `Protocol deviation` = 100 * mean(pd),
                                `BBS` = 100 * mean(bbs), `Genetic disease` = 100 * mean(genetic)) %>%
  pivot_longer(-c(route, n)) %>% mutate(name = factor(name, levels = c("Protocol deviation", "Retinal dystrophy", "Genetic disease", "BBS")))
pS2 <- ggplot(rt, aes(name, value, fill = sprintf("%s (n=%d)", route, n))) +
  geom_col(position = position_dodge(.8), width = .8, colour = "black", linewidth = .2) +
  geom_text(aes(label = sprintf("%.1f", value)), position = position_dodge(.8), vjust = -.4, size = 3) +
  scale_fill_manual(values = c("black", "grey75")) + scale_y_continuous(limits = c(0, 108), expand = c(0, 0)) +
  labs(x = NULL, y = "% of adults (≥18 years)") + theme_w
ggsave(file.path(FIG, "FigureS2_adults_by_route.png"), pS2, width = 7.5, height = 4.5, dpi = 300)
message("06 figures done")
