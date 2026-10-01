# Purpose: Primary Cox models (SAP Sec. 6.3) and secondary outcomes (SAP Sec. 4; review B3).
# Input(s): data/derived/analysis.rds
# Output(s): output/tables/table2_allele_models.csv     (Model 1-3; per-allele HR; 2-df joint test)
#            output/tables/s08_genotype_models.csv      (6 genotypes vs e3/e3 and E4 carrier; Model 3)
#            output/tables/s08_ph_test.csv              (Schoenfeld residual tests)
#            output/tables/s08_time_split.csv           (only if PH violated for an APOE term)
#            output/tables/s08_secondary_outcomes.csv   (Model 3 allele model)
#            output/models/primary_m3.rds
# Key methods/packages: survival::coxph (Efron ties), survival::cox.zph, survival::survSplit
# Notes: Primary test (proposed, review A1) = Model 3 2-df LRT of e2 and e4 counts.
#        With common_sample = true, Models 1-3 share the Model 3 complete-case sample.
#        kinship_method = "cluster" switches to robust SE and a robust Wald joint test.

source("R/00_setup.R")
log_msg("== 08 primary Cox models ==")
a <- load_derived("analysis", cfg)
cluster <- if (identical(cfg$cohort$kinship_method, "cluster")) "family_id" else NULL
if (!is.null(cluster)) a[is.na(family_id), family_id := id]
common <- if (isTRUE(cfg$models$common_sample)) model_covs(cfg, "m3") else NULL

# ---- Table 2: allele-additive models, Model 1-3 -----------------------------------
fits <- lapply(c("m1", "m2", "m3"), function(m)
  fit_apoe(cfg, a, "primary", m, cluster = cluster, sample_vars = common,
           label = c(m1 = "Model 1: age, sex", m2 = "Model 2: + CKD risk factors",
                     m3 = "Model 3: + lipids")[[m]]))
t2 <- data.table::rbindlist(lapply(fits, `[[`, "table"))
t2[, `:=`(term_label = term_label(term), `HR (95% CI)` = fmt_ci(hr, lcl, ucl), p_fmt = fmt_p(p),
          joint_p_fmt = fmt_p(joint_p))]
write_table(t2, "table2_allele_models", cfg)
saveRDS(fits[[3]]$fit, file.path(cfg$paths$output, "models", "primary_m3.rds"))

# ---- Genotype categories and E4 carrier, Model 3 -------------------------------------
gm <- fit_apoe(cfg, a, "primary", "m3", exposure = "apoe_genotype", cluster = cluster,
               sample_vars = common, label = "Genotype (ref e3/e3), Model 3")$table
ec <- fit_apoe(cfg, a, "primary", "m3", exposure = "e4_carrier", cluster = cluster,
               sample_vars = common, label = "E4 carrier vs non-E4, Model 3")$table
gtab <- rbind(gm, ec)
gtab[, `:=`(term_label = term_label(term), `HR (95% CI)` = fmt_ci(hr, lcl, ucl), p_fmt = fmt_p(p))]
write_table(gtab, "s08_genotype_models", cfg)

# ---- Proportional hazards check (Model 3) ----------------------------------------------
zph <- survival::cox.zph(fits[[3]]$fit, transform = "km")
ph <- data.table::data.table(term = rownames(zph$table), chisq = zph$table[, "chisq"],
                             df = zph$table[, "df"], p = zph$table[, "p"])
write_table(ph, "s08_ph_test", cfg)
viol <- ph[term %in% c("e2_count", "e4_count") & p < cfg$models$ph_alpha, term]

# ---- Time-split estimates if PH is violated for an APOE term -------------------------------
if (length(viol)) {
  log_msg("PH violated for: ", paste(viol, collapse = ", "), " -> time-split model")
  dat <- fits[[3]]$data
  cuts <- unlist(cfg$models$time_split_years)
  sp <- survival::survSplit(Surv(tstart, primary_tstop, primary_event) ~ ., data = as.data.frame(dat),
                            cut = cuts, episode = "period")
  sp$period <- factor(sp$period, labels = c(paste0("0-", cuts[1]),
                                            paste0(cuts[-length(cuts)], "-", cuts[-1]),
                                            paste0(">", cuts[length(cuts)])))
  covs <- model_covs(cfg, "m3")
  f <- stats::as.formula(paste("Surv(tstart, primary_tstop, primary_event) ~ e2_count:period + e4_count:period +",
                               paste(covs, collapse = " + "), if (!is.null(cluster)) "+ cluster(family_id)"))
  ft <- survival::coxph(f, data = sp, ties = "efron")
  terms <- grep("^e[24]_count:period", names(stats::coef(ft)), value = TRUE)
  ts <- tidy_hr(ft, terms)
  ts[, `HR (95% CI)` := fmt_ci(hr, lcl, ucl)]
} else {
  ts <- data.table::data.table(note = "Proportional hazards not rejected for e2_count or e4_count; time-split model not fitted.")
}
write_table(ts, "s08_time_split", cfg)

# ---- Secondary outcomes, Model 3 allele model -------------------------------------------------
sec <- setdiff(outcome_names(), "primary")
st <- data.table::rbindlist(lapply(sec, function(o) {
  dat <- if (o == "egfr_lt60") a[egfr_index >= cfg$outcome$egfr_threshold] else a
  res <- tryCatch(fit_apoe(cfg, dat, o, "m3", cluster = cluster, sample_vars = common,
                           label = outcome_labels()[[o]])$table,
                  error = function(e) data.table::data.table(label = outcome_labels()[[o]], outcome = o,
                                                             note = conditionMessage(e)))
  res
}), fill = TRUE)
if ("hr" %in% names(st)) st[!is.na(hr), `:=`(term_label = term_label(term), `HR (95% CI)` = fmt_ci(hr, lcl, ucl))]
write_table(st, "s08_secondary_outcomes", cfg)

m3 <- t2[model == "m3"]
log_msg("08 done: Model 3 n = ", m3$n[1], ", events = ", m3$events[1], ", e2 HR ",
        m3[term == "e2_count", `HR (95% CI)`], ", e4 HR ", m3[term == "e4_count", `HR (95% CI)`],
        ", 2-df P = ", m3$joint_p_fmt[1])
