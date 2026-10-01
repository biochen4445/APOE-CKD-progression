# Input data schema

Place one CSV per table in `data/raw/` (real EMR extracts, never committed) or let `sim/simulate_data.R` write synthetic versions to `data/simulated/`. Column names are case-sensitive. Dates are `YYYY-MM-DD`. `id` is the study pseudo-ID shared by all tables.

## person.csv (one row per participant)

| Column | Required | Description |
| --- | --- | --- |
| id | yes | Study ID |
| sex | yes | `M` or `F` |
| birth_date | yes | Date of birth |
| last_visit_date | yes | Last CMUH encounter of any type |
| death_date | yes (may be empty) | Date of death; source to be confirmed (in-hospital only, or linked national death registry) |
| consent_date | optional | iHi biobank consent date; needed for delayed-entry analysis (review A3) |
| family_id | optional | Family / pedigree ID; needed for kinship analyses (review A2) |
| death_cause_code, death_cause_system | optional | Underlying cause of death code and its system (`ICD9CM` / `ICD10CM`) |

## apoe_genotype.csv (one row per genotyped participant)

| Column | Required | Description |
| --- | --- | --- |
| id | yes | Study ID |
| rs429358 | yes | Count or dosage (0-2) of the e4-defining allele (C) |
| rs7412 | yes | Count or dosage (0-2) of the e2-defining allele (T) |
| info_rs429358, info_rs7412 | optional | Imputation INFO score |
| pc1 ... pc10 | optional | Genetic principal components |

## labs.csv (long format)

| Column | Description |
| --- | --- |
| id, date | Study ID, specimen date |
| test | One of `creatinine` (mg/dL), `hdl`, `ldl`, `tg` (mg/dL), `hba1c` (%), `glucose_fasting`, `glucose_random` (mg/dL) |
| value | Numeric result in the units above |
| setting | `OPD`, `IPD` or `ER` |

## diagnoses.csv

| Column | Description |
| --- | --- |
| id, date | Study ID; encounter date (admission date for inpatient stays) |
| encounter_id | Admission ID (required for `IPD` rows; may be empty otherwise) |
| discharge_date | Discharge date for `IPD` rows (optional; admission date used if missing) |
| setting | `OPD`, `IPD` or `ER` |
| code_system | `ICD9CM` or `ICD10CM` |
| code | Diagnosis code, with or without dots |
| dx_rank | 1 = principal diagnosis, 2+ = secondary |

## procedures.csv

| Column | Description |
| --- | --- |
| id, date, setting | As above |
| encounter_id | Admission ID for inpatient procedures |
| code_system | `ICD9PROC`, `ICD10PCS` or `NHI` |
| code | Procedure / order code |

## medications.csv

| Column | Description |
| --- | --- |
| id, date | Study ID, prescription date |
| atc | ATC code |
| days_supply | optional |

## vitals.csv

| Column | Description |
| --- | --- |
| id, date, setting | As above |
| sbp, dbp | mm Hg |
| height_cm, weight_kg | May be empty |

## catastrophic.csv (optional)

Catastrophic illness certificates (重大傷病), if linked: `id`, `date` (certificate start), `category` (`ESRD` used for ESKD).
