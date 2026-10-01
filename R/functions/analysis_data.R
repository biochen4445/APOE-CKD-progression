# Purpose: Merge cohort, APOE, outcomes and covariates into one analysis data set, and build
# follow-up time variables. `make_analysis_data()` reruns stages 03-05 with a modified config,
# which the sensitivity analyses (10_sensitivity.R) use.

assemble_analysis <- function(cfg, outcomes, covs, apoe) {
  keep_out <- c("id", "index_date", "censor_date", "died", "death_date",
                intersect(c("consent_date", "family_id"), names(outcomes)),
                grep(paste0("^(", paste(outcome_names(), collapse = "|"), ")_(date|event|exit|status)$"),
                     names(outcomes), value = TRUE),
                "cr_rise_date", "ckd_hosp_date", "ckd_death_date")
  keep_out <- unique(intersect(keep_out, names(outcomes)))
  a <- merge(covs, outcomes[, ..keep_out], by = "id")
  pcs <- grep("^pc[0-9]+$", names(apoe), value = TRUE)
  a <- merge(a, apoe[, c("id", "apoe_genotype", "e2_count", "e4_count", pcs), with = FALSE], by = "id")
  a[, apoe_genotype := factor(apoe_genotype, levels = unlist(cfg$apoe$genotype_order))]
  a[, apoe_genotype := stats::relevel(apoe_genotype, ref = cfg$apoe$reference_genotype)]
  a[, e4_carrier := as.integer(e4_count > 0)]

  # Entry time: index date, or max(index, consent) with delayed entry (review A3)
  a[, entry_date := index_date]
  if (isTRUE(cfg$cohort$delayed_entry)) {
    stop_if(!"consent_date" %in% names(a), "delayed_entry = true requires person.consent_date")
    a[!is.na(consent_date), entry_date := as_idate(pmax(index_date, consent_date))]
  }
  a[, tstart := as.numeric(entry_date - index_date) / 365.25]
  for (o in outcome_names()) {
    a[, (paste0(o, "_tstop")) := as.numeric(get(paste0(o, "_exit")) - index_date) / 365.25]
    a[, (paste0(o, "_age_exit")) := age + get(paste0(o, "_tstop"))]
    # with delayed entry, records whose exit precedes entry are not at risk for that outcome
    a[get(paste0(o, "_tstop")) <= tstart,
      c(paste0(o, "_tstop"), paste0(o, "_event"), paste0(o, "_status")) := list(NA_real_, NA_integer_, NA_integer_)]
  }
  a[, age_entry := age + tstart]
  a[, fu_years := primary_tstop - tstart]
  a[]
}

make_analysis_data <- function(cfg, d, apoe, cl) {
  coh <- build_cohort(cfg, d, apoe, cl)
  out <- derive_outcomes(cfg, d, coh$cohort, cl)
  cov <- derive_covariates(cfg, d, coh$cohort, cl)
  list(analysis = assemble_analysis(cfg, out, cov, apoe), flow = coh$flow)
}

# Deep-modify a config list: modify_cfg(cfg, list(cohort = list(index_setting = "ANY")))
modify_cfg <- function(cfg, changes) utils::modifyList(cfg, changes)
