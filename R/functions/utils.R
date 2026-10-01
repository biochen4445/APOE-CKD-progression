# Purpose: General helpers shared by all stage scripts (logging, I/O, code matching, windows).
# Key packages: data.table

log_msg <- function(...) {
  msg <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(...))
  message(msg)
  lf <- getOption("apoeckd.logfile")
  if (!is.null(lf)) cat(msg, "\n", file = lf, append = TRUE, sep = "")
  invisible(msg)
}

stop_if <- function(cond, ...) if (isTRUE(cond)) stop(paste0(...), call. = FALSE)

# Restore integer storage of IDate after pmin/pmax/ifelse (which can return double storage)
as_idate <- function(x) structure(as.integer(round(unclass(x))), class = c("IDate", "Date"))

# Age in completed fractional years
age_years <- function(birth_date, date) as.numeric(date - birth_date) / 365.25

# Normalise diagnosis / procedure / ATC codes: strip dots and spaces, upper case
norm_code <- function(x) toupper(gsub("[^A-Za-z0-9]", "", x))

# --- I/O -----------------------------------------------------------------------
raw_dir <- function(cfg) {
  if (identical(cfg$data_source, "emr")) cfg$paths$raw_emr else cfg$paths$raw_simulated
}

read_raw <- function(cfg, name, required = TRUE) {
  f <- file.path(raw_dir(cfg), paste0(name, ".csv"))
  if (!file.exists(f)) {
    stop_if(required, "Required input not found: ", f)
    return(NULL)
  }
  data.table::fread(f, colClasses = list(character = "id"), na.strings = c("", "NA"))
}

save_derived <- function(x, name, cfg) {
  saveRDS(x, file.path(cfg$paths$derived, paste0(name, ".rds")))
}
load_derived <- function(name, cfg) {
  f <- file.path(cfg$paths$derived, paste0(name, ".rds"))
  stop_if(!file.exists(f), "Derived file missing: ", f, " -- run the earlier stage first.")
  readRDS(f)
}

write_table <- function(x, name, cfg) {
  f <- file.path(cfg$paths$output, "tables", paste0(name, ".csv"))
  data.table::fwrite(x, f)
  log_msg("Wrote ", f)
  invisible(f)
}

save_figure <- function(p, name, cfg, width = 7, height = 4.5) {
  f <- file.path(cfg$paths$output, "figures", paste0(name, ".pdf"))
  ggplot2::ggsave(f, p, width = width, height = height)
  ggplot2::ggsave(sub("\\.pdf$", ".png", f), p, width = width, height = height, dpi = 300)
  log_msg("Wrote ", f)
  invisible(f)
}

# --- Code lists -------------------------------------------------------------------
load_codelists <- function(cfg) {
  cl <- data.table::fread(file.path(cfg$paths$codelists, "codelists.csv"),
                          colClasses = "character")
  cl[, code_prefix := norm_code(code_prefix)]
  stop_if(any(cl$code_prefix == ""), "Empty code_prefix in codelists.csv would match every code.")
  cl
}

# TRUE where `code` (in `system`) starts with any prefix of list `list_name`
match_codelist <- function(code, system, list_name, codelists) {
  out <- logical(length(code))
  sub <- codelists[list == list_name]
  stop_if(nrow(sub) == 0, "Unknown code list: ", list_name)
  for (s in unique(sub$code_system)) {
    pat <- paste0("^(", paste(sub[code_system == s, code_prefix], collapse = "|"), ")")
    idx <- which(system == s)
    out[idx] <- grepl(pat, code[idx])
  }
  out
}

# --- Windows around the index date --------------------------------------------------
# Closest record to index within [index + lo, index + hi] days; ties go to the earlier record.
closest_in_window <- function(events, anchor, value_cols, lo, hi, date_col = "date") {
  x <- merge(events, anchor[, .(id, index_date)], by = "id")
  x[, d := as.numeric(get(date_col) - index_date)]
  x <- x[d >= lo & d <= hi]
  x[, absd := abs(d)]
  data.table::setorderv(x, c("id", "absd", "d"))
  x <- x[!duplicated(id)]
  x[, c("id", value_cols), with = FALSE]
}

# Any record within [index + lo, index + hi] days
any_in_window <- function(events, anchor, lo, hi, date_col = "date") {
  x <- merge(events, anchor[, .(id, index_date)], by = "id")
  x[, d := as.numeric(get(date_col) - index_date)]
  unique(x[d >= lo & d <= hi, id])
}

# Exact Poisson CI for a rate (events / person-time), Garwood
pois_exact_ci <- function(events, pt, level = 0.95) {
  a <- (1 - level) / 2
  lo <- ifelse(events == 0, 0, stats::qchisq(a, 2 * events) / 2)
  hi <- stats::qchisq(1 - a, 2 * (events + 1)) / 2
  data.table::data.table(rate = events / pt, lcl = lo / pt, ucl = hi / pt)
}

fmt_ci <- function(est, lo, hi, digits = 2) {
  sprintf(paste0("%.", digits, "f (%.", digits, "f-%.", digits, "f)"), est, lo, hi)
}
fmt_p <- function(p) ifelse(is.na(p), NA_character_,
                            ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
