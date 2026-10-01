# Purpose: Load packages, configuration and helper functions; create output folders.
# Input(s): config/config.yml, config/codelists/codelists.csv, R/functions/*.R
# Output(s): objects `cfg` and `cl` in the calling environment; output/ and data/derived/ folders
# Key methods/packages: data.table, survival, ggplot2, yaml
# Notes: Every stage script sources this file, so each stage can also be run on its own.
#        Run from the repository root (the folder containing config/).

if (!exists(".apoeckd_setup_done") || !isTRUE(.apoeckd_setup_done)) {
  if (!file.exists("config/config.yml")) {
    stop("Working directory must be the repository root (folder containing config/).", call. = FALSE)
  }
  suppressPackageStartupMessages({
    library(data.table)
    library(survival)
    library(ggplot2)
    library(yaml)
  })

  cfg <- yaml::read_yaml("config/config.yml", fileEncoding = "UTF-8")
  for (f in sort(list.files("R/functions", pattern = "\\.R$", full.names = TRUE))) source(f)

  for (p in c(cfg$paths$derived, file.path(cfg$paths$output, c("tables", "figures", "logs", "models")))) {
    dir.create(p, recursive = TRUE, showWarnings = FALSE)
  }
  if (is.null(getOption("apoeckd.logfile"))) {
    options(apoeckd.logfile = file.path(cfg$paths$output, "logs",
                                        paste0("pipeline_", format(Sys.Date(), "%Y%m%d"), ".log")))
  }
  set.seed(cfg$seed)
  cl <- load_codelists(cfg)
  .apoeckd_setup_done <- TRUE
  log_msg("Setup complete. data_source = ", cfg$data_source, "; SAP ", cfg$sap_version)
}
