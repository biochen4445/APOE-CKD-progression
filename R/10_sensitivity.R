# Purpose: Sensitivity analyses (SAP Sec. 6.5) plus review-proposed analyses (docs/SAP_review_v0.1_HL.md).
# Input(s): data/derived/{inputs,apoe,analysis}.rds
# Output(s): output/tables/s10_sensitivity.csv; output/figures/s10_sensitivity_forest.pdf/.png
# Key methods/packages: Model 3 allele-additive Cox model refitted under each scenario
# Notes: Scenarios that change the cohort or outcome definition rebuild stages 03-05 with a
#        modified config (make_analysis_data). Review scenarios run only when the needed
#        columns exist (consent_date, family_id, pc1-pc10).

source("R/00_setup.R")
log_msg("== 10 sensitivity analyses ==")
d <- load_derived("inputs", cfg)
apoe <- load_derived("apoe", cfg)
a <- load_derived("analysis", cfg)

run <- function(label, source, dat, cfg_use = cfg, ...) {
  log_msg("  scenario: ", label)
  res <- tryCatch(fit_apoe(cfg_use, dat, "primary", "m3", label = label, ...)$table,
                  error = function(e) data.table::data.table(label = label, note = conditionMessage(e)))
  res[, source := source]
  res
}
rebuild <- function(changes) make_analysis_data(modify_cfg(cfg, changes), d, apoe, cl)$analysis

sc <- list()
sc[[1]]  <- run("Primary analysis (Model 3)", "Primary", a)
sc[[2]]  <- run("Age at index 45-64 y (ARIC range)", "SAP 6.5", a[age < 65])
sc[[3]]  <- run("Index creatinine from any setting", "SAP 6.5", rebuild(list(cohort = list(index_setting = "ANY"))))
sc[[4]]  <- run("MDRD eGFR as covariate", "SAP 6.5", data.table::copy(a)[, egfr_index := egfr_index_mdrd])
sc[[5]]  <- run("Creatinine rise without 90-day confirmation", "SAP 6.5",
                rebuild(list(outcome = list(require_confirmation = FALSE))))
sc[[6]]  <- run("Exclude e2/e4 genotype", "SAP 6.5", a[apoe_genotype != "e2/e4"])
sc[[7]]  <- run("Age as time scale", "SAP 6.5", a, timescale = "age")
sc[[8]]  <- run("Exclude index eGFR < 60", "SAP 6.5", a[egfr_index >= 60])
sc[[9]]  <- run("Exclude prevalent CKD diagnosis code", "Review A6", a[ckd_dx_prevalent == 0])
sc[[10]] <- run("Confirmation by any later value", "Review A7",
                rebuild(list(outcome = list(confirm_rule = "any"))))
if ("consent_date" %in% names(a) && any(!is.na(a$consent_date))) {
  sc[[11]] <- run("Delayed entry at biobank consent", "Review A3",
                  rebuild(list(cohort = list(delayed_entry = TRUE))))
}
if ("family_id" %in% names(a) && any(!is.na(a$family_id))) {
  sc[[12]] <- run("One participant per family", "Review A2",
                  rebuild(list(cohort = list(kinship_method = "unrelated"))))
  ac <- data.table::copy(a)[is.na(family_id), family_id := id]
  sc[[13]] <- run("Robust SE clustered on family", "Review A2", ac, cluster = "family_id")
}
if (all(paste0("pc", seq_len(cfg$models$n_pcs)) %in% names(a))) {
  sc[[14]] <- run("Adjusted for genetic PC1-10", "Review C", a,
                  cfg_use = modify_cfg(cfg, list(models = list(adjust_pcs = TRUE))))
}

sens <- data.table::rbindlist(sc, fill = TRUE)
if ("hr" %in% names(sens)) sens[!is.na(hr), `:=`(term_label = term_label(term), `HR (95% CI)` = fmt_ci(hr, lcl, ucl),
                                                 joint_p_fmt = fmt_p(joint_p))]
write_table(sens, "s10_sensitivity", cfg)

pd <- sens[!is.na(hr)]
pd[, label := factor(label, levels = rev(unique(label)))]
pd[, source_group := factor(sub(" .*", "", source), levels = c("Primary", "SAP", "Review"),
                            labels = c("Primary", "SAP 6.5", "Review-proposed"))]
p <- ggplot(pd, aes(x = hr, y = label)) +
  geom_vline(xintercept = 1, colour = "grey60") +
  geom_errorbarh(aes(xmin = lcl, xmax = ucl), height = 0.2, colour = "grey30") +
  geom_point(aes(shape = source_group), size = 2.2, colour = "#1f3b6f") +
  facet_wrap(~ term_label, nrow = 1) +
  scale_x_log10(breaks = c(0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.25, 1.5, 2)) +
  labs(x = "Hazard ratio per allele (95% CI), Model 3, log scale", y = NULL, shape = "Source") +
  theme_minimal(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom")
save_figure(p, "s10_sensitivity_forest", cfg, width = 9, height = 5.5)
log_msg("10 done: ", length(unique(sens$label)), " scenarios")
