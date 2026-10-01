# Purpose: Index date, inclusion and exclusion criteria (SAP Sec. 2), with a step-by-step flow table.
# Index date = first serum creatinine at age [age_min, age_max] (inclusive), outpatient by default.

build_cohort <- function(cfg, d, apoe, cl) {
  cc <- cfg$cohort
  flow <- data.table::data.table(step = character(), n = integer(), note = character())
  add_flow <- function(step, ids, note = "") {
    flow <<- rbind(flow, data.table::data.table(step = step, n = length(unique(ids)), note = note))
  }

  person <- d$person
  add_flow("iHi participants with EMR person record", person$id)

  # Inclusion 1: APOE genotype passing QC
  ids <- intersect(person$id, apoe[apoe_qc_pass == TRUE, id])
  add_flow("I1: APOE genotype available and passing QC", ids)

  # Exclusion 3 (applied here because age needs birth date): missing age or sex
  p <- person[id %in% ids]
  p <- p[!is.na(birth_date) & sex %in% c("M", "F")]
  add_flow("E3: excluded missing birth date (age) or sex", p$id, "SAP lists as exclusion 3; applied before index search")

  # Inclusion 2: >= 1 creatinine at age 45-65
  cr <- d$labs[test == "creatinine" & !is.na(value) & value > 0 & id %in% p$id]
  if (toupper(cc$index_setting) != "ANY") cr <- cr[setting == cc$index_setting]
  cr <- merge(cr, p[, .(id, birth_date)], by = "id")
  cr[, age := age_years(birth_date, date)]
  cr <- cr[age >= cc$age_min & age < cc$age_max + 1]
  idx <- cr[, .(index_cr = mean(value), index_age = age[1], index_setting = setting[1]),
            by = .(id, date)]
  data.table::setorder(idx, id, date)
  idx <- idx[!duplicated(id)]
  data.table::setnames(idx, "date", "index_date")
  add_flow(sprintf("I2: >= 1 %s creatinine at age %d-%d", cc$index_setting, cc$age_min, cc$age_max),
           idx$id, "same-day values averaged")

  # Inclusion 3: >= 1 encounter after index (followable)
  idx <- merge(idx, p[, .(id, sex, birth_date, last_visit_date)], by = "id")
  idx <- idx[last_visit_date > index_date]
  add_flow("I3: >= 1 encounter after index date", idx$id)

  # Exclusion 1: severe kidney dysfunction at index (ARIC thresholds)
  idx <- idx[!((sex == "M" & index_cr >= cc$severe_cr_male) |
               (sex == "F" & index_cr >= cc$severe_cr_female))]
  add_flow(sprintf("E1: excluded index creatinine >= %.1f (M) / %.1f (F) mg/dL",
                   cc$severe_cr_male, cc$severe_cr_female), idx$id)

  # Exclusion 2: dialysis, kidney transplant or ESKD on or before index
  dx <- d$diagnoses[id %in% idx$id]
  dx <- dx[match_codelist(code, code_system, "eskd_baseline_dx", cl)]
  pr <- d$procedures[id %in% idx$id]
  pr <- pr[match_codelist(code, code_system, "dialysis_proc", cl) |
           match_codelist(code, code_system, "transplant_proc", cl)]
  prior <- unique(c(
    any_in_window(dx, idx, -Inf, 0),
    any_in_window(pr, idx, -Inf, 0)
  ))
  if (!is.null(d$catastrophic)) {
    prior <- unique(c(prior, any_in_window(d$catastrophic[category == "ESRD"], idx, -Inf, 0)))
  }
  idx <- idx[!id %in% prior]
  add_flow("E2: excluded dialysis, kidney transplant or ESKD on or before index", idx$id)

  # Optional (not in SAP v0.1-HL): keep one member per family
  if (identical(cc$kinship_method, "unrelated")) {
    fam <- person[id %in% idx$id, .(id, family_id)]
    idx <- merge(idx, fam, by = "id")
    idx[is.na(family_id), family_id := id]
    data.table::setorder(idx, family_id, index_date, id)
    idx <- idx[!duplicated(family_id)]
    idx[, family_id := NULL]
    add_flow("Optional: one participant per family retained", idx$id, "review A2; not in SAP v0.1-HL")
  }

  add_flow("Analytic cohort", idx$id)
  list(cohort = idx[], flow = flow[])
}
