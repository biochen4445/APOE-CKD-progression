# Purpose: Derive primary and secondary outcomes, censoring and follow-up (SAP Sec. 4, Sec. 5).
# Input(s): data/derived/inputs.rds, data/derived/cohort.rds
# Output(s): data/derived/outcomes.rds; output/tables/s04_outcome_counts.csv,
#            s04_primary_components.csv
# Key methods/packages: data.table rolling joins for confirmed (sustained) laboratory events
# Notes: Creatinine rise requires confirmation by the next OPD value >= 90 days later (config).
#        CKD hospitalisation event date = discharge date (JAMA 2005), admission if missing.

source("R/00_setup.R")
log_msg("== 04 outcomes ==")
d <- load_derived("inputs", cfg)
cohort <- load_derived("cohort", cfg)

out <- derive_outcomes(cfg, d, cohort, cl)
save_derived(out, "outcomes", cfg)

lab <- outcome_labels()
counts <- data.table::rbindlist(lapply(outcome_names(), function(o) {
  at_risk <- if (o == "egfr_lt60") out[egfr_index_primary >= cfg$outcome$egfr_threshold] else out
  data.table::data.table(outcome = o, label = lab[[o]], n_at_risk = nrow(at_risk),
                         events = sum(at_risk[[paste0(o, "_event")]]),
                         deaths_without_event = sum(at_risk[[paste0(o, "_status")]] == 2L))
}))

# Which component defined the primary event date (JAMA 2005 reported this split)
comp <- out[primary_event == 1, .(
  cr_first   = !is.na(cr_rise_date) & cr_rise_date == primary_date,
  hosp_first = !is.na(ckd_hosp_date) & ckd_hosp_date == primary_date,
  death_first = !is.na(ckd_death_date) & ckd_death_date == primary_date,
  cr_ever = !is.na(cr_rise_date), hd_ever = !is.na(ckd_hosp_date) | !is.na(ckd_death_date))]
components <- data.table::data.table(
  component = c("Creatinine rise only", "CKD hospitalisation/death only", "Both criteria",
                "Event date set by creatinine rise", "Event date set by hospitalisation",
                "Event date set by CKD death"),
  n = c(comp[cr_ever & !hd_ever, .N], comp[!cr_ever & hd_ever, .N], comp[cr_ever & hd_ever, .N],
        comp[cr_first == TRUE, .N], comp[hosp_first == TRUE, .N], comp[death_first == TRUE, .N]))

write_table(counts, "s04_outcome_counts", cfg)
write_table(components, "s04_primary_components", cfg)
log_msg("04 done: primary events = ", sum(out$primary_event))
