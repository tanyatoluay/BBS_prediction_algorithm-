# =============================================================================
# 01 — Build the patient-level analysis dataset (N = 1004)
# =============================================================================

# ---- 1. eCRF V4.0 -------------------------------------------------------------
# Two header rows; several labels are duplicated (checklist vs organ-system section),
# so columns are addressed by spreadsheet position and given explicit names.
ecrf_raw <- read_excel(ECRF_FILE, sheet = "Final Report", col_names = FALSE, skip = 2,
                       col_types = "text", .name_repair = "minimal")
stopifnot(ncol(ecrf_raw) == 75)
col_map <- c(
  pd = 1, site = 2, id = 3, sex = 8, ethnicity = 9, height = 11, weight = 12,
  # eligibility checklist (paediatric section)
  ck_rodcone = 16, ck_renal = 17, ck_ataxia = 18, ck_syn = 19, ck_poly = 20, ck_brachy = 21,
  ck_hyperphagia = 22, ck_cog = 23, ck_speech = 24, ck_hypog = 25,
  # eligibility checklist (adult section)
  ad_bmi30 = 26, ad_obese_lt6 = 27, ad_retinal = 28, clin_dx_bbs = 29,
  # organ-system section (Ja / Nein / Unklar)
  os_visual = 30, os_retdeg = 31, os_optic = 32, os_conerod = 33, os_rp = 34,
  os_renalfail = 35, os_fanconi = 36, os_tubulo = 37, os_fsgs = 38, os_nephrotic = 39, os_protein = 40,
  os_hypog = 41, os_cardiomyo = 42, os_hypert = 43, os_hyperpig = 44,
  os_poly = 45, os_syn = 46, os_brachy = 47, os_crowding = 48, os_hypodontia = 49, os_micropenis = 50,
  os_cog = 51, os_seizure = 52, os_myoclon = 53, os_motor = 54, os_decline = 55, os_speech = 56,
  os_ataxia = 57, os_id = 58,
  age_months = 63, lt192 = 64,
  # derived columns in the export (BN–BR are headed "ERN criteria", BS–BW "Beales criteria")
  exp_ern_prim_n = 66, exp_ern_sec_n = 68, exp_ern_met = 70,
  exp_beales_prim_n = 71, exp_beales_prim_sym = 72, exp_beales_sec_n = 73, exp_beales_met = 75
)
ecrf <- ecrf_raw[, col_map]; names(ecrf) <- names(col_map)
ecrf <- ecrf %>% mutate(across(everything(), ~ str_squish(.x)))

norm_yn <- function(x) {
  x <- tolower(str_squish(x))
  case_when(x %in% c("ja", "yes") ~ "yes", x %in% c("nein", "no") ~ "no",
            x %in% c("unklar", "unclear") ~ "unclear", TRUE ~ NA_character_)
}
ecrf <- ecrf %>% mutate(across(c(pd, starts_with("ck_"), starts_with("ad_"), starts_with("os_"),
                                 clin_dx_bbs), norm_yn))

# Three-state feature: present (any item yes) / not assessed / absent.
#   not_assessed rule "unclear": >=1 item recorded "unclear" and none "yes"  (manuscript definition)
tri <- function(df, cols) {
  m <- as.matrix(df[, cols, drop = FALSE])
  any_yes <- rowSums(m == "yes", na.rm = TRUE) > 0
  any_unc <- rowSums(m == "unclear", na.rm = TRUE) > 0
  any_no  <- rowSums(m == "no", na.rm = TRUE) > 0
  factor(case_when(any_yes ~ "present", any_unc ~ "not_assessed",
                   any_no ~ "absent", TRUE ~ "not_assessed"),
         levels = c("present", "absent", "not_assessed"))
}
FEATURE_ITEMS <- list(
  retinal     = c("ck_rodcone", "ad_retinal", "os_conerod", "os_retdeg", "os_rp"),
  visual      = c("os_visual", "os_optic"),
  poly        = c("ck_poly", "os_poly"),
  syn         = c("ck_syn", "os_syn"),
  brachy      = c("ck_brachy", "os_brachy"),
  renal       = c("ck_renal", "os_renalfail", "os_fanconi", "os_tubulo", "os_fsgs", "os_nephrotic", "os_protein"),
  hypog       = c("ck_hypog", "os_hypog"),
  micropenis  = c("os_micropenis"),
  cog         = c("ck_cog", "os_cog"),
  id_dis      = c("os_id"),
  decline     = c("os_decline"),
  speech      = c("ck_speech", "os_speech"),
  motor       = c("os_motor"),
  ataxia      = c("ck_ataxia", "os_ataxia"),
  seizures    = c("os_seizure", "os_myoclon"),
  hyperphagia = c("ck_hyperphagia"),
  dental      = c("os_crowding", "os_hypodontia"),
  hypert      = c("os_hypert"),
  cardiomyo   = c("os_cardiomyo"),
  hyperpig    = c("os_hyperpig")
)
FEATURE_ITEMS$hypog_genital <- c(FEATURE_ITEMS$hypog, FEATURE_ITEMS$micropenis)
FEATURE_ITEMS$synbrachy     <- c(FEATURE_ITEMS$syn, FEATURE_ITEMS$brachy)
FEATURE_ITEMS$nd            <- c(FEATURE_ITEMS$cog, FEATURE_ITEMS$id_dis, FEATURE_ITEMS$decline,
                                 FEATURE_ITEMS$speech, FEATURE_ITEMS$motor)
FEATURE_ITEMS$ophth         <- c(FEATURE_ITEMS$retinal, FEATURE_ITEMS$visual)

for (f in names(FEATURE_ITEMS)) {
  ecrf[[paste0(f, "_3")]] <- tri(ecrf, FEATURE_ITEMS[[f]])
  ecrf[[f]] <- as.integer(ecrf[[paste0(f, "_3")]] == "present")   # primary: not assessed = absent
}

# ---- 2. Genetic database V21.0 ---------------------------------------------------
grp_raw <- read_excel(GEN_FILE, sheet = "groups")
grp_lab <- c("BBS", "Non-BBS obesity-associated", "Non-obesity", "No qualifying diagnosis")
groups <- map2_dfr(seq_along(grp_raw), grp_lab, ~ tibble(id = str_squish(as.character(grp_raw[[.x]])), group = .y)) %>%
  filter(!is.na(id), id != "")
stopifnot(!anyDuplicated(groups$id), nrow(groups) == 1004)

gen_all <- read_excel(GEN_FILE, sheet = "all data", col_types = "text") %>%
  mutate(cg_number = str_squish(cg_number))
# header labels contain non-breaking spaces -> locate by pattern
names(gen_all)[str_detect(names(gen_all), "BMI-SDS")]         <- "eCRF – BMI-SDS (AGA/KHS) Kinder"
names(gen_all)[str_detect(names(gen_all), "BMI Erwachsene")]  <- "eCRF – BMI Erwachsene (kg/m²)"
names(gen_all)[str_detect(names(gen_all), "^Reclassification primary")] <-
  "Reclassification primary 29.09.2026 ARO (V18.6, automated + manual)"
per_id <- gen_all %>% filter(cg_number %in% groups$id) %>%
  group_by(cg_number) %>%
  summarise(
    relation = paste(sort(unique(relation_to_index)), collapse = "|"),
    bmi_sds  = suppressWarnings(as.numeric(first(na.omit(`eCRF – BMI-SDS (AGA/KHS) Kinder`)))),
    bmi_adult = suppressWarnings(as.numeric(first(na.omit(`eCRF – BMI Erwachsene (kg/m²)`)))),
    .groups = "drop") %>%
  mutate(index = as.integer(relation == "index"))   # 983 index; 21 relatives fulfilling eligibility

# ---- 3. Merge -----------------------------------------------------------------------
dat <- groups %>% inner_join(ecrf, by = "id") %>% left_join(per_id, by = c("id" = "cg_number"))
stopifnot(nrow(dat) == 1004)

dat <- dat %>% mutate(
  group      = factor(group, levels = grp_lab),
  grp4       = as.integer(group),
  bbs        = as.integer(group == "BBS"),
  genetic    = as.integer(grp4 %in% 1:3),
  obesity_gd = as.integer(grp4 %in% 1:2),
  site       = str_pad(site, 2, pad = "0"),
  female     = as.integer(sex %in% c("weiblich", "Weiblich")),
  male       = 1L - female,
  white      = as.integer(ethnicity == "Weiß"),
  pd         = as.integer(pd %in% "yes"),
  clin_dx    = as.integer(clin_dx_bbs %in% "yes"),
  age_months = suppressWarnings(as.numeric(age_months)),
  age_y      = age_months / 12,
  age18      = factor(if_else(age_months < 216, "<18", ">=18"), levels = c("<18", ">=18")),
  age16      = factor(if_else(age_months < 192, "<16", ">=16"), levels = c("<16", ">=16")),
  bbs_cat    = case_when(id %in% BBS_CAT$A ~ "A", id %in% BBS_CAT$B1 ~ "B1",
                         id %in% BBS_CAT$B2 ~ "B2", id %in% BBS_CAT$C ~ "C", TRUE ~ NA_character_),
  high_site  = as.integer(site %in% HIGH_BBS_SITES),
  # BMI from eCRF height/weight (reproduces V25 table 1; V21 adult-BMI column is incomplete)
  bmi_calc   = suppressWarnings(as.numeric(weight) / (as.numeric(height) / 100)^2),
  # eligibility route actually documented in the eCRF
  route      = case_when(!is.na(ad_retinal) & is.na(ck_poly) ~ "adult section only",
                         !is.na(ad_retinal) & !is.na(ck_poly) ~ "both sections",
                         !is.na(ck_poly) ~ "paediatric checklist only", TRUE ~ "none"),
  hyperphagia_assessed = !is.na(ck_hyperphagia)
)
# data-quality flags (to be queried with PBE / sites)
dq_flags <- dat %>% filter(is.na(age_months) | age_months < 0) %>% select(id, site, group, age_months)
write_xlsx(dq_flags, file.path(OUT_DIR, "tables", "check_data_quality_age.xlsx"))
stopifnot(sum(!is.na(dat$bbs_cat)) == 38, all(dat$bbs[!is.na(dat$bbs_cat)] == 1))

# ---- 4. Composite variables, criteria, Score A --------------------------------------
dat <- dat %>% mutate(
  learning = as.integer(cog | id_dis | decline),
  # domains (methods: seven predefined organ-system domains)
  dom_eye   = as.integer(retinal | visual),
  dom_cns   = as.integer(nd | ataxia | seizures),
  dom_renal = renal,
  dom_endo  = as.integer(hypog | micropenis | hyperpig),
  dom_app   = hyperphagia,
  dom_skel  = as.integer(poly | syn | brachy | dental),
  dom_cv    = as.integer(cardiomyo | hypert),
  n_domains = dom_eye + dom_cns + dom_renal + dom_endo + dom_app + dom_skel + dom_cv,
  dom3      = as.integer(n_domains >= 3),
  n_features = retinal + visual + poly + syn + brachy + renal + hypog + micropenis + cog + id_dis +
               decline + speech + motor + ataxia + seizures + hyperphagia + dental + hypert + cardiomyo + hyperpig,
  # Beales (all eCRF items; obesity counted for everybody by eligibility)
  beales_prim = 1 + retinal + poly + learning + hypog_genital + renal,
  beales_sec  = speech + synbrachy + motor + ataxia + dental + cardiomyo,
  beales      = as.integer(beales_prim >= 4 | (beales_prim == 3 & beales_sec >= 2)),
  # V25 manuscript variant: syndactyly and brachydactyly counted as TWO secondary features
  # (reproduces 89 / 63 / 47 in V25; contradicts the V25 methods and Beales 1999)
  beales_sec_v25 = speech + syn + brachy + motor + ataxia + dental + cardiomyo,
  beales_v25  = as.integer(beales_prim >= 4 | (beales_prim == 3 & beales_sec_v25 >= 2)),
  beales_count = beales_prim + beales_sec / 10,       # ordinal count, secondary as tie-breaker
  # ERN (clinical items only)
  ern_prim = 1 + poly + retinal + renal,
  ern_sec  = hypog_genital + nd,
  ern      = as.integer(ern_prim >= 4 | (ern_prim == 3 & ern_sec >= 2)),
  # export-derived criteria (eligibility checklist only)
  beales_export = as.integer(exp_beales_met == "Yes"),
  ern_export    = as.integer(exp_ern_met == "Yes"),
  # Score A
  scoreA = 5 * poly + 2 * retinal + 2 * hypog_genital + 2 * renal + synbrachy + nd,
  scoreA_noret = scoreA - 2 * retinal,
  scoreA3 = as.integer(scoreA >= 3), scoreA5 = as.integer(scoreA >= 5)
)

# Score A with not-assessed items counted as present
pres_or_na <- function(f) as.integer(dat[[paste0(f, "_3")]] != "absent")
dat$scoreA_napres <- 5 * pres_or_na("poly") + 2 * pres_or_na("retinal") + 2 * pres_or_na("hypog_genital") +
  2 * pres_or_na("renal") + pres_or_na("synbrachy") + pres_or_na("nd")
dat$scoreA_complete <- with(dat, poly_3 != "not_assessed" & retinal_3 != "not_assessed" &
  hypog_genital_3 != "not_assessed" & renal_3 != "not_assessed" & synbrachy_3 != "not_assessed" &
  nd_3 != "not_assessed")

# ---- 5. Cross-check of BBS evidence categories against V21.0 VUS tiers -------------------
bbs_vus <- gen_all %>% filter(cg_number %in% unlist(BBS_CAT), str_detect(finding_primary_class, "Class 3")) %>%
  transmute(id = cg_number, gene = finding_primary_gene, zyg = finding_primary_zygosity,
            tier_V21_colJ = str_extract(vus_description_primary, "VUS-(high|mid|low)"),
            tier_V18_6_ARO = str_extract(`Reclassification primary 29.09.2026 ARO (V18.6, automated + manual)`,
                                         "^(P|LP|VUS-high|VUS-mid|VUS-low|VUS|LB|B)\\b")) %>%
  left_join(dat %>% select(id, bbs_cat, beales, ern), by = "id")
write_xlsx(bbs_vus, file.path(OUT_DIR, "tables", "check_BBS_VUS_tiers.xlsx"))

saveRDS(dat, file.path(OUT_DIR, "wegio_analysis_dataset.rds"))
message("Dataset built: ", nrow(dat), " participants; BBS ", sum(dat$bbs))
