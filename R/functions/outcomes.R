# Purpose: Outcome derivation (SAP Sec. 4) and follow-up / censoring (SAP Sec. 5).
# Primary: CKD progression = earliest of (1) confirmed creatinine rise >= 0.4 mg/dL from index,
# (2) CKD hospitalisation (CKD code any position, principal diagnosis not AKI), (3) death with CKD.
# Censoring = earliest of death, last CMUH encounter, EMR end date.

# First date at which `flag` is TRUE and confirmed >= confirm_days later.
# x: data.table(id, date, flag) restricted to the at-risk window.
# rule "next": the first measurement >= confirm_days after the candidate must also be flagged.
# rule "any" : any measurement >= confirm_days after the candidate is flagged.
# A candidate with no measurement >= confirm_days later is not an event (review A7).
detect_sustained <- function(x, confirm_days = 90, rule = c("next", "any"),
                             require_confirmation = TRUE) {
  rule <- match.arg(rule)
  x <- data.table::copy(x)[!is.na(flag)]
  data.table::setorder(x, id, date)
  if (!require_confirmation) return(x[flag == TRUE, .(event_date = min(date)), by = id])
  if (rule == "any") {
    x[, conf := rev(cummax(rev(as.integer(flag)))) == 1L, by = id]
  } else {
    x[, conf := flag]
  }
  cand <- x[flag == TRUE, .(id, cand_date = date, target = date + as.integer(confirm_days))]
  if (nrow(cand) == 0) return(data.table::data.table(id = character(), event_date = x$date[0]))
  meas <- x[, .(id, target = date, conf)]
  data.table::setkey(meas, id, target)
  j <- meas[cand, on = .(id, target), roll = -Inf]
  j[conf == TRUE, .(event_date = min(cand_date)), by = id]
}

# Attach an event date column to `co` (named <name>_date), keeping only dates in (start, censor]
attach_event <- function(co, ev, name) {
  col <- paste0(name, "_date")
  co[, (col) := as_idate(ev$event_date[match(co$id, ev$id)])]
  co[!is.na(get(col)) & (get(col) <= index_date | get(col) > censor_date), (col) := NA]
  co
}

derive_outcomes <- function(cfg, d, cohort, cl) {
  oc <- cfg$outcome
  co <- data.table::copy(cohort)
  co <- merge(co, d$person[, intersect(c("id", "death_date", "consent_date", "family_id",
                                         "death_cause_code", "death_cause_system"),
                                       names(d$person)), with = FALSE], by = "id")
  emr_end <- data.table::as.IDate(cfg$emr_end_date)
  co[, censor_date := as_idate(pmin(death_date, last_visit_date, emr_end, na.rm = TRUE))]
  co[, died := !is.na(death_date) & death_date <= censor_date]

  # --- (1) Creatinine rise --------------------------------------------------------
  cr <- d$labs[test == "creatinine" & !is.na(value) & id %in% co$id]
  if (toupper(oc$cr_rise_setting) != "ANY") cr <- cr[setting == oc$cr_rise_setting]
  cr <- merge(cr[, .(id, date, value)], co[, .(id, index_date, index_cr, censor_date, sex, birth_date)], by = "id")
  cr <- cr[date > index_date & date <= censor_date]
  cr[, flag := value - index_cr >= oc$cr_rise_mgdl]
  ev <- detect_sustained(cr[, .(id, date, flag)], oc$confirm_days, oc$confirm_rule,
                         isTRUE(oc$require_confirmation))
  co <- attach_event(co, ev, "cr_rise")

  # --- (2) CKD hospitalisation -------------------------------------------------------
  dx <- d$diagnoses[setting == "IPD" & id %in% co$id]
  dx[, ckd := match_codelist(code, code_system, "ckd_outcome_dx", cl)]
  dx[, aki_principal := dx_rank == 1 & match_codelist(code, code_system, "aki_primary_exclude_dx", cl)]
  enc <- dx[, .(ckd = any(ckd), aki_principal = any(aki_principal),
                admit = min(date),
                discharge = suppressWarnings(max(discharge_date, na.rm = TRUE))),
            by = .(id, encounter_id)]
  enc[!is.finite(discharge), discharge := NA]
  pr <- d$procedures[setting == "IPD" & id %in% co$id & !is.na(encounter_id)]
  pr_ckd <- unique(pr[match_codelist(code, code_system, "ckd_outcome_proc", cl), .(id, encounter_id)])
  enc[pr_ckd, ckd := TRUE, on = .(id, encounter_id)]
  enc[, event_date := if (identical(oc$hosp_event_date, "discharge"))
    data.table::fifelse(is.na(discharge), admit, discharge) else admit]
  hosp <- merge(enc[ckd & !aki_principal], co[, .(id, index_date)], by = "id")[admit > index_date]
  ev <- hosp[, .(event_date = min(event_date)), by = id]
  co <- attach_event(co, ev, "ckd_hosp")

  # --- (3) Death with CKD -----------------------------------------------------------
  dd <- co[died == TRUE, .(id, death_date)]
  ckd_any <- d$diagnoses[id %in% dd$id]
  ckd_any <- ckd_any[match_codelist(code, code_system, "ckd_outcome_dx", cl), .(id, date)]
  ckd_any <- merge(ckd_any, dd, by = "id")
  death_ckd <- unique(ckd_any[date <= death_date & date >= death_date - oc$death_dx_window_days, id])
  if ("death_cause_code" %in% names(co)) {
    dc <- co[died == TRUE & !is.na(death_cause_code)]
    death_ckd <- unique(c(death_ckd,
      dc[match_codelist(norm_code(death_cause_code), death_cause_system, "ckd_outcome_dx", cl), id]))
  }
  co[, ckd_death_date := as_idate(data.table::fifelse(id %in% death_ckd, death_date, data.table::as.IDate(NA)))]
  co[ckd_death_date <= index_date, ckd_death_date := NA]

  # --- Primary composite and its components -------------------------------------------
  co[, primary_date := as_idate(pmin(cr_rise_date, ckd_hosp_date, ckd_death_date, na.rm = TRUE))]
  co[, ckd_hosp_death_date := as_idate(pmin(ckd_hosp_date, ckd_death_date, na.rm = TRUE))]

  # --- Secondary eGFR-based outcomes (OPD creatinine, primary equation) ------------
  cr[, egfr := egfr_calc(value, age_years(birth_date, date), sex, cfg$egfr$primary_equation,
                         cfg$egfr$mdrd_constant)]
  co[, egfr_index_primary := egfr_calc(index_cr, index_age, sex, cfg$egfr$primary_equation,
                                       cfg$egfr$mdrd_constant)]
  cr <- merge(cr, co[, .(id, egfr_index_primary)], by = "id")
  cr[, flag := egfr < oc$egfr_threshold]
  ev <- detect_sustained(cr[egfr_index_primary >= oc$egfr_threshold, .(id, date, flag)],
                         oc$confirm_days, oc$confirm_rule, TRUE)
  co <- attach_event(co, ev, "egfr_lt60")
  cr[, flag := egfr <= (1 - oc$egfr_decline_fraction) * egfr_index_primary]
  ev <- detect_sustained(cr[, .(id, date, flag)], oc$confirm_days, oc$confirm_rule, TRUE)
  co <- attach_event(co, ev, "egfr_decline40")

  # --- ESKD: chronic dialysis (record >= 90 d after first) or transplant ----------------
  dia <- rbind(
    d$procedures[id %in% co$id & match_codelist(code, code_system, "dialysis_proc", cl), .(id, date)],
    d$diagnoses[id %in% co$id & match_codelist(code, code_system, "dialysis_dx", cl), .(id, date)])
  dia <- merge(unique(dia), co[, .(id, index_date, censor_date)], by = "id")[date > index_date & date <= censor_date]
  dia[, flag := TRUE]
  ev_dia <- detect_sustained(dia[, .(id, date, flag)], oc$eskd_dialysis_days, "any", TRUE)
  tx <- rbind(
    d$procedures[id %in% co$id & match_codelist(code, code_system, "transplant_proc", cl), .(id, date)],
    d$diagnoses[id %in% co$id & match_codelist(code, code_system, "transplant_dx", cl), .(id, date)])
  tx <- merge(tx, co[, .(id, index_date)], by = "id")[date > index_date]
  ev_tx <- first_event(tx)
  ev_eskd <- rbind(ev_dia, ev_tx)
  if (!is.null(d$catastrophic)) {
    ci <- merge(d$catastrophic[category == "ESRD", .(id, date)], co[, .(id, index_date)], by = "id")
    ev_eskd <- rbind(ev_eskd, first_event(ci[date > index_date]))
  }
  ev_eskd <- ev_eskd[, .(event_date = min(event_date)), by = id]
  co <- attach_event(co, ev_eskd, "eskd")

  # --- Event indicators and follow-up time -----------------------------------------------
  for (o in outcome_names()) {
    dcol <- paste0(o, "_date")
    co[, (paste0(o, "_event")) := as.integer(!is.na(get(dcol)))]
    co[, (paste0(o, "_exit")) := as_idate(data.table::fifelse(is.na(get(dcol)), censor_date, get(dcol)))]
    # competing-risk status: 0 censored, 1 event, 2 death without the event
    co[, (paste0(o, "_status")) := data.table::fifelse(!is.na(get(dcol)), 1L,
                                    data.table::fifelse(died, 2L, 0L))]
  }
  co[]
}

# Earliest date per id (empty-safe)
first_event <- function(x) {
  if (!nrow(x)) return(data.table::data.table(id = character(), event_date = data.table::as.IDate(integer())))
  x[, .(event_date = min(date)), by = id]
}

outcome_names <- function() c("primary", "cr_rise", "ckd_hosp_death", "egfr_lt60", "eskd", "egfr_decline40")

outcome_labels <- function() c(
  primary        = "CKD progression (composite, primary)",
  cr_rise        = "Creatinine rise >= 0.4 mg/dL",
  ckd_hosp_death = "CKD hospitalisation or death",
  egfr_lt60      = "Incident eGFR < 60 (index eGFR >= 60)",
  eskd           = "ESKD (chronic dialysis or transplant)",
  egfr_decline40 = "Sustained eGFR decline >= 40%")
