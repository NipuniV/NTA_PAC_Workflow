# =============================================================================
# 01_helpers.R — Functions shared by more than one workflow
# -----------------------------------------------------------------------------
# Keeping these in one place guarantees that, for example, WF1 and WF3 reshape
# spectra using identical logic, and that WF2 and the KMD stage compute Kendrick
# mass defects with identical constants. In the original separate projects each
# script carried its own copy of these functions, which is how small differences
# (such as the aromaticity-index fold threshold) were able to creep in.
# =============================================================================


# -----------------------------------------------------------------------------
# ppm_diff() — parts-per-million difference between an observed and a reference
# mass. The reference is used as the denominator, which is the convention when
# comparing against a tabulated (known) value.
# -----------------------------------------------------------------------------
ppm_diff <- function(obs, ref) abs(obs - ref) / ref * 1e6


# -----------------------------------------------------------------------------
# formula_to_mono() — monoisotopic exact mass from a molecular formula string.
# Handles multi-letter symbols (Cl, Br, Si) and implicit counts ("CH4" = C1 H4).
# Returns NA if the formula contains an element not in ATOMIC_MONO.
#
# Used by 02_validate_libraries.R to check every reference table against itself.
# -----------------------------------------------------------------------------
formula_to_mono <- function(formula) {
  if (is.na(formula) || formula == "") return(NA_real_)
  # each token = an element symbol (capital + optional lowercase) + optional count
  m <- str_match_all(formula, "([A-Z][a-z]?)(\\d*)")[[1]]
  if (nrow(m) == 0) return(NA_real_)
  els <- m[, 2]
  cnt <- ifelse(m[, 3] == "", 1L, suppressWarnings(as.integer(m[, 3])))
  if (any(!els %in% names(ATOMIC_MONO))) return(NA_real_)
  sum(ATOMIC_MONO[els] * cnt)
}


# -----------------------------------------------------------------------------
# read_spectra() — load a ChromaTOF wide-format export.
#
# Expected column layout:
#   1  Name
#   2  1st Dimension Time (s)
#   3  2nd Dimension Time (s)
#   4  Base Mass
#   5+ alternating mass / relative-intensity pairs: m1, ri1, m2, ri2, ...
#      with intensities normalised so the base peak = 1000
#
# clean_names() converts these to: name, x1st_dimension_time_s,
# x2nd_dimension_time_s, base_mass, m1, ri1, m2, ri2, ...
# -----------------------------------------------------------------------------
read_spectra <- function(path) {
  stopifnot(file.exists(path))
  clean_names(read_csv(path, show_col_types = FALSE))
}


# -----------------------------------------------------------------------------
# reshape_ms_long() — pivot the wide m/ri columns into one row per fragment ion.
#
# id_cols identifies a unique peak. Relative intensity is recomputed as a
# percentage of the largest ion in that spectrum so that spectra acquired at
# different absolute intensities are directly comparable.
#
# top_n = NULL keeps every ion (needed by WF3, which tests all ion pairs for
# neutral losses). top_n = 4 keeps only the four most intense ions (WF1's
# spectral fingerprint check).
# -----------------------------------------------------------------------------
reshape_ms_long <- function(df, id_cols, top_n = NULL) {
  out <- df %>%
    pivot_longer(
      cols          = matches("^(m|ri)\\d+$"),
      names_to      = c(".value", "peak_num"),
      names_pattern = "([a-z]+)(\\d+)"
    ) %>%
    filter(!is.na(m), m > 0, !is.na(ri), ri > 0) %>%
    group_by(across(all_of(id_cols))) %>%
    mutate(rel_int_pct = ri / max(ri) * 100) %>%
    ungroup()

  if (!is.null(top_n)) {
    out <- out %>%
      group_by(across(all_of(id_cols))) %>%
      slice_max(rel_int_pct, n = top_n, with_ties = FALSE) %>%
      arrange(desc(rel_int_pct), .by_group = TRUE) %>%
      mutate(ion_rank = row_number()) %>%
      ungroup()
  }
  out
}


# -----------------------------------------------------------------------------
# remove_column_bleed() — drop peaks whose BASE MASS matches a known
# polysiloxane column-bleed ion within `ppm_tol`.
#
# Bleed peaks are chromatographically ubiquitous and would otherwise be carried
# through every filter as spurious candidates.
# -----------------------------------------------------------------------------
remove_column_bleed <- function(df, ppm_tol) {
  df %>%
    rowwise() %>%
    mutate(.is_bleed = any(ppm_diff(base_mass, bleed_masses) <= ppm_tol)) %>%
    ungroup() %>%
    filter(!.is_bleed) %>%
    select(-.is_bleed)
}


# -----------------------------------------------------------------------------
# aromaticity_index() — AI computed directly from a measured ion mass.
#
#   1. take the decimal part of the mass (the mass defect)
#   2. fold it to a symmetric scale: if the decimal exceeds `fold`, subtract 1
#      so that values just below the next integer become small negatives
#   3. normalise by mass/13 (13 = nominal mass of a CH unit)
#
# Aromatic, hydrogen-poor ions cluster in a narrow band close to zero;
# hydrogen-rich aliphatics sit well above it.
#
# The fold threshold and the acceptance window are set in 00_setup.R.
# -----------------------------------------------------------------------------
aromaticity_index <- function(mass, fold = wf2_ai_fold) {
  dec    <- mass - floor(mass)
  defect <- ifelse(dec > fold, dec - 1, dec)
  defect / (mass / 13)
}


# -----------------------------------------------------------------------------
# add_kmd() — add Kendrick Mass and Kendrick Mass Defect columns for every
# series defined in `kmd_series` (CH2, H/Cl, Br2).
#
# Creates columns km_CH2 / kmd_CH2, km_Cl / kmd_Cl, km_Br2 / kmd_Br2 from the
# column named by `mass_col`.
# -----------------------------------------------------------------------------
add_kmd <- function(df, mass_col = "exact_mass") {
  for (ser in names(kmd_series)) {
    e <- kmd_series[[ser]]$exact
    n <- kmd_series[[ser]]$nominal
    km  <- df[[mass_col]] * (n / e)
    df[[paste0("km_",  ser)]] <- km
    df[[paste0("kmd_", ser)]] <- round(km) - km
  }
  df
}


# -----------------------------------------------------------------------------
# extract_mono_peak() — pull the monoisotopic (lowest) mass out of a formatted
# isotope-peak string such as:
#     "167.97863 (ri=1000), 169.97557 (ri=367)"
# The minimum value is the M peak; the higher values are M+2 / M+4
# isotopologues. Using M is the chemically correct basis for KMD.
# -----------------------------------------------------------------------------
extract_mono_peak <- function(peak_str) {
  if (is.na(peak_str)) return(NA_real_)
  vals <- suppressWarnings(
    as.numeric(str_extract_all(peak_str, "[0-9]+\\.[0-9]+")[[1]])
  )
  if (length(vals) == 0 || all(is.na(vals))) NA_real_ else min(vals, na.rm = TRUE)
}


# -----------------------------------------------------------------------------
# stage_file() — copy a produced file into the data folder of the workflow that
# consumes it next, so each data_wfN folder always contains everything that
# workflow needs to run on its own.
# -----------------------------------------------------------------------------
stage_file <- function(from, to_dir) {
  if (file.exists(from)) {
    file.copy(from, file.path(to_dir, basename(from)), overwrite = TRUE)
    invisible(TRUE)
  } else {
    warning("stage_file: missing ", from)
    invisible(FALSE)
  }
}

message("Helpers loaded.")
