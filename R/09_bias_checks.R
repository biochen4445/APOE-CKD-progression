# Purpose: Bias checks (SAP Sec. 6.4): follow-up and mortality by genotype; death as a competing risk
#          (cause-specific vs Fine-Gray); creatinine testing intensity by genotype.
# Input(s): data/derived/analysis.rds, data/derived/inputs.rds
# Output(s): output/tables/s09_followup_death_by_genotype.csv, s09_competing_risk.csv,
#            s09_testing_intensity.csv
# Key methods/packages: survival::finegray + coxph (subdistribution HR); quasi-Poisson regression
# Notes: The cause-specific HR equals the primary Cox model (deaths censored).

source("R/00_setup.R")
log_msg("== 09 bias checks ==")
a <- load_derived("analysis", cfg)
d <- load_derived("inputs", cfg)
a <- a[!is.na(primary_tstop)]
a[, py := primary_tstop - tstart]
glev <- unlist(cfg$apoe$genotype_order)
common <- if (isTRUE(cfg$models$common_sample)) model_covs(cfg, "m3") else NULL

# ---- 1. Follow-up and death by genotype ---------------------------------------------------
fd <- data.table::rbindlist(lapply(c(glev, "Overall"), function(g) {
  x <- if (g == "Overall") a else a[apoe_genotype == g]
  ci <- pois_exact_ci(sum(x$died), sum(x$py))
  data.table::data.table(genotype = g, n = nrow(x),
                         median_fu_y = round(stats::median(x$py), 2),
                         iqr_fu = sprintf("%.1f-%.1f", stats::quantile(x$py, .25), stats::quantile(x$py, .75)),
                         deaths = sum(x$died),
                         death_rate_per_1000py = round(1000 * ci$rate, 2),
                         death_rate_ci = sprintf("%.2f-%.2f", 1000 * ci$lcl, 1000 * ci$ucl))
}))
dm <- fit_apoe(cfg, data.table::copy(a)[, death_event := as.integer(died)][, death_tstop := primary_tstop],
               "death", "m1", label = "All-cause death, Model 1")$table
fd_model <- dm[, .(genotype = term_label(term), death_HR_model1 = fmt_ci(hr, lcl, ucl), p = fmt_p(p))]
write_table(rbind(fd, fd_model, fill = TRUE), "s09_followup_death_by_genotype", cfg)

# ---- 2. Competing risk: cause-specific vs Fine-Gray (Model 3) -------------------------------
covs <- model_covs(cfg, "m3")
dat <- complete_sample(a, unique(c("tstart", "primary_tstop", "primary_status", "e2_count", "e4_count", covs, common)))
csh <- fit_apoe(cfg, dat, "primary", "m3", label = "Cause-specific HR (deaths censored)")$table
dat[, status_f := factor(primary_status, levels = 0:2, labels = c("censor", "ckd", "death"))]
delayed <- any(dat$tstart > 0)
fg_lhs <- if (delayed) "Surv(tstart, primary_tstop, status_f)" else "Surv(primary_tstop, status_f)"
fgd <- survival::finegray(stats::as.formula(paste(fg_lhs, "~ .")),
                          data = as.data.frame(dat[, c(if (delayed) "tstart", "primary_tstop", "status_f",
                                                       "e2_count", "e4_count", covs), with = FALSE]),
                          etype = "ckd")
fgf <- survival::coxph(stats::as.formula(paste("Surv(fgstart, fgstop, fgstatus) ~",
                                              paste(c("e2_count", "e4_count", covs), collapse = " + "))),
                       data = fgd, weights = fgwt)
sh <- tidy_hr(fgf, c("e2_count", "e4_count"))
sh[, `:=`(label = "Fine-Gray subdistribution HR", outcome = "primary", model = "m3",
          n = nrow(dat), events = sum(dat$primary_status == 1))]
cr_tab <- rbind(csh[, .(label, term, n, events, hr, lcl, ucl, p)], sh[, .(label, term, n, events, hr, lcl, ucl, p)])
cr_tab[, `:=`(term_label = term_label(term), `HR (95% CI)` = fmt_ci(hr, lcl, ucl), deaths_competing = sum(dat$primary_status == 2))]
write_table(cr_tab, "s09_competing_risk", cfg)

# ---- 3. Creatinine testing intensity during follow-up ---------------------------------------
cr <- d$labs[test == "creatinine" & id %in% a$id]
if (toupper(cfg$outcome$cr_rise_setting) != "ANY") cr <- cr[setting == cfg$outcome$cr_rise_setting]
a[, exit_date := index_date + as.integer(round(primary_tstop * 365.25))]
cnt <- merge(cr[, .(id, date)], a[, .(id, entry_date, exit_date)], by = "id")[date > entry_date & date <= exit_date, .N, by = id]
a[, n_tests := cnt$N[match(a$id, cnt$id)]]
a[is.na(n_tests), n_tests := 0L]
ti <- a[py > 0, .(n = .N, tests_per_py_mean = round(sum(n_tests) / sum(py), 3),
                  median_tests = stats::median(n_tests)), by = apoe_genotype][order(match(apoe_genotype, glev))]
qp <- stats::glm(stats::as.formula(paste("n_tests ~ e2_count + e4_count +", paste(model_covs(cfg, "m1"), collapse = " + "))),
                 family = stats::quasipoisson(), offset = log(py), data = a[py > 0])
b <- stats::coef(summary(qp))[c("e2_count", "e4_count"), ]
rr <- data.table::data.table(apoe_genotype = c("Rate ratio per e2 allele (Model 1)", "Rate ratio per e4 allele (Model 1)"),
                             tests_per_py_mean = sprintf("%.3f (%.3f-%.3f)", exp(b[, 1]),
                                                         exp(b[, 1] - 1.96 * b[, 2]), exp(b[, 1] + 1.96 * b[, 2])),
                             median_tests = fmt_p(b[, 4]))
data.table::setnames(rr, "median_tests", "p")
write_table(rbind(ti[, .(apoe_genotype = as.character(apoe_genotype), n, tests_per_py_mean = as.character(tests_per_py_mean),
                         median_tests)], rr, fill = TRUE), "s09_testing_intensity", cfg)
log_msg("09 done")
