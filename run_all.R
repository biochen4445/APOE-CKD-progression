# Purpose: Run the full pipeline (stages 01-11) and write an audit log with key numbers,
#          configuration fingerprint, git commit and package versions.
# Usage:   Rscript run_all.R            (from the repository root)
#          Rscript run_all.R --simulate (regenerate synthetic data first)
# Output(s): output/tables/*, output/figures/*, output/logs/run_log.md, output/logs/sessionInfo.txt

args <- commandArgs(trailingOnly = TRUE)
if ("--simulate" %in% args) source("sim/simulate_data.R")

source("R/00_setup.R")
t0 <- Sys.time()
stages <- sort(list.files("R", pattern = "^[0-9]{2}_.*\\.R$", full.names = TRUE))
stages <- stages[!grepl("^R/00_", stages)]
timing <- data.table::data.table(stage = basename(stages), seconds = NA_real_)
for (k in seq_along(stages)) {
  started <- Sys.time()
  source(stages[k], local = new.env(parent = globalenv()))  # each stage in its own environment
  timing$seconds[k] <- round(as.numeric(difftime(Sys.time(), started, units = "secs")), 1)
}

# ---- Audit log ------------------------------------------------------------------------
tab <- function(n) data.table::fread(file.path(cfg$paths$output, "tables", paste0(n, ".csv")))
flow <- tab("fig0_cohort_flow"); t2 <- tab("table2_allele_models"); oc <- tab("s04_outcome_counts")
git <- tryCatch(system("git rev-parse --short HEAD", intern = TRUE, ignore.stderr = TRUE),
                error = function(e) NA, warning = function(w) NA)
cfg_md5 <- unname(tools::md5sum("config/config.yml"))
cl_md5 <- unname(tools::md5sum("config/codelists/codelists.csv"))
writeLines(capture.output(sessionInfo()), file.path(cfg$paths$output, "logs", "sessionInfo.txt"))

pk <- c("data.table", "survival", "ggplot2", "yaml")
pkv <- paste(sprintf("%s %s", pk, vapply(pk, function(p) as.character(utils::packageVersion(p)), "")), collapse = "; ")
m3 <- t2[model == "m3"]
log <- c(
  "# Run log -- APOE x CKD progression",
  "",
  sprintf("- Run: %s (%.1f min)", format(t0, "%Y-%m-%d %H:%M"), as.numeric(difftime(Sys.time(), t0, units = "mins"))),
  sprintf("- Data source: **%s**%s", cfg$data_source,
          if (cfg$data_source == "simulated") " (synthetic data; numbers are NOT study results)" else ""),
  sprintf("- SAP: %s", cfg$sap_version),
  sprintf("- Git commit: %s", if (length(git) && !is.na(git[1])) git[1] else "not a git repository / uncommitted"),
  sprintf("- config.yml md5: %s; codelists.csv md5: %s", cfg_md5, cl_md5),
  sprintf("- %s; %s", R.version.string, pkv),
  "",
  "## Key numbers (pulled from output/tables)",
  "",
  sprintf("- Analytic cohort: %s", format(flow[step == "Analytic cohort", n], big.mark = ",")),
  sprintf("- Primary events: %s", oc[outcome == "primary", events]),
  sprintf("- Model 3 sample: n = %s, events = %s", m3$n[1], m3$events[1]),
  sprintf("- Model 3 per e2 allele HR: %s; per e4 allele HR: %s; 2-df %s P = %s",
          m3[term == "e2_count", `HR (95% CI)`], m3[term == "e4_count", `HR (95% CI)`],
          m3$joint_type[1], m3$joint_p_fmt[1]),
  "",
  "## Stage timing (seconds)",
  "",
  paste0("- ", timing$stage, ": ", timing$seconds),
  "",
  "## Open items",
  "",
  "- See docs/SAP_review_v0.1_HL.md; code lists marked `draft-verify` need checking before the real run."
)
writeLines(log, file.path(cfg$paths$output, "logs", "run_log.md"))
log_msg("Pipeline complete. See output/logs/run_log.md")
