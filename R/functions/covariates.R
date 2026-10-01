# Purpose: Index covariates (SAP Sec. 5). Windows are days relative to index (config$windows).
#   Labs: closest value in [-365, +30]; diagnoses: any time on or before index;
#   medications: any prescription in [-90, 0]; BP and BMI: closest OPD record in [-180, +180].

derive_covariates <- function(cfg, d, cohort, cl) {
  w <- cfg$windows
  co <- data.table::copy(cohort[, .(id, index_date, index_age, index_cr, sex)])
  out <- co[, .(id, age = index_age, sex = factor(sex, levels = c("F", "M")))]

  # eGFR at index: primary equation and MDRD (sensitivity Sec. 6.5)
  out[, egfr_index := egfr_calc(co$index_cr, co$index_age, co$sex, cfg$egfr$primary_equation,
                                cfg$egfr$mdrd_constant)]
  out[, egfr_index_mdrd := egfr_mdrd(co$index_cr, co$index_age, co$sex, cfg$egfr$mdrd_constant)]
  out[, index_cr := co$index_cr]

  # Lipids: closest value in the lab window
  labs <- d$labs[id %in% co$id & !is.na(value)]
  for (t in c("hdl", "ldl", "tg")) {
    v <- closest_in_window(labs[test == t, .(id, date, value)], co, "value", -w$lab_before, w$lab_after)
    out[, (t) := v$value[match(out$id, v$id)]]
  }
  out[, log_tg := log(tg)]

  # Diabetes: diagnosis on/before index, OR any qualifying lab in the lab window, OR A10 in the drug window
  dm_dx <- any_in_window(d$diagnoses[id %in% co$id & match_codelist(code, code_system, "diabetes_dx", cl)],
                         co, -Inf, 0)
  lab_win <- function(test_name, thr) {
    any_in_window(labs[test == test_name & value >= thr], co, -w$lab_before, w$lab_after)
  }
  dm_lab <- unique(c(lab_win("hba1c", w$dm_hba1c), lab_win("glucose_fasting", w$dm_fpg),
                     lab_win("glucose_random", w$dm_random_glucose)))
  meds <- d$medications[id %in% co$id]
  dm_rx <- any_in_window(meds[match_codelist(atc, rep("ATC", length(atc)), "glucose_lowering_atc", cl)],
                         co, -w$med_before, 0)
  out[, diabetes := as.integer(id %in% c(dm_dx, dm_lab, dm_rx))]

  # Antihypertensive use in the drug window
  htn_rx <- any_in_window(meds[match_codelist(atc, rep("ATC", length(atc)), "antihtn_atc", cl)],
                          co, -w$med_before, 0)
  out[, antihtn := as.integer(id %in% htn_rx)]

  # Coronary heart disease history: diagnosis or revascularisation on/before index
  chd <- unique(c(
    any_in_window(d$diagnoses[id %in% co$id & match_codelist(code, code_system, "chd_dx", cl)], co, -Inf, 0),
    any_in_window(d$procedures[id %in% co$id & match_codelist(code, code_system, "revasc_proc", cl)], co, -Inf, 0)))
  out[, chd := as.integer(id %in% chd)]

  # Blood pressure and BMI: closest OPD record with the measure present
  vit <- d$vitals[id %in% co$id & setting == "OPD"]
  bp <- closest_in_window(vit[!is.na(sbp) & !is.na(dbp), .(id, date, sbp, dbp)], co,
                          c("sbp", "dbp"), -w$vitals_halfwidth, w$vitals_halfwidth)
  out[, sbp := bp$sbp[match(out$id, bp$id)]]
  out[, dbp := bp$dbp[match(out$id, bp$id)]]
  bm <- closest_in_window(vit[!is.na(height_cm) & !is.na(weight_kg), .(id, date, height_cm, weight_kg)],
                          co, c("height_cm", "weight_kg"), -w$vitals_halfwidth, w$vitals_halfwidth)
  bm[, bmi := weight_kg / (height_cm / 100)^2]
  bm[bmi < 12 | bmi > 70, bmi := NA]
  out[, bmi := bm$bmi[match(out$id, bm$id)]]

  # Prevalent CKD diagnosis code on/before index (review A6 sensitivity)
  ckd_prev <- any_in_window(d$diagnoses[id %in% co$id & match_codelist(code, code_system, "ckd_outcome_dx", cl)],
                            co, -Inf, 0)
  out[, ckd_dx_prevalent := as.integer(id %in% ckd_prev)]

  out[]
}
