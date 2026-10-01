# Purpose: Define index date and apply inclusion / exclusion criteria (SAP Sec. 2); Figure 0 data.
# Input(s): data/derived/inputs.rds, data/derived/apoe.rds
# Output(s): data/derived/cohort.rds; output/tables/fig0_cohort_flow.csv, fig0_cohort_flow.md
# Key methods/packages: data.table
# Notes: Index = first creatinine at age 45-65 (inclusive) in the configured setting (OPD primary).

source("R/00_setup.R")
log_msg("== 03 cohort ==")
d <- load_derived("inputs", cfg)
apoe <- load_derived("apoe", cfg)

res <- build_cohort(cfg, d, apoe, cl)
res$flow[, excluded := c(NA, -diff(n))]
save_derived(res$cohort, "cohort", cfg)
write_table(res$flow, "fig0_cohort_flow", cfg)

# Mermaid flow diagram (renders on GitHub)
fl <- res$flow
lines <- c("```mermaid", "flowchart TD")
for (i in seq_len(nrow(fl))) {
  lines <- c(lines, sprintf('  S%d["%s<br/>n = %s"]', i, gsub('"', "'", fl$step[i]),
                            format(fl$n[i], big.mark = ",")))
  if (i > 1) lines <- c(lines, sprintf("  S%d --> S%d", i - 1, i))
}
lines <- c(lines, "```")
writeLines(lines, file.path(cfg$paths$output, "tables", "fig0_cohort_flow.md"))
log_msg("03 done: analytic cohort n = ", nrow(res$cohort))
