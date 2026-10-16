# =============================================================================
# 02_validate_libraries.R — Automated reference-library integrity check
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
#
# Every workflow depends on tabulated reference values. If a tabulated exact
# mass is wrong, the affected compound simply fails to match and disappears
# silently from the results — there is no error message, and the loss is very
# easy to miss.
#
# This actually happened: the reference entry for 1,6-dibromopyrene carried
# 360.0419 Da, which is the AVERAGE molecular weight of C16H8Br2, not its
# monoisotopic mass (357.8993 Da). The isotope filter correctly recognised the
# compound as di-brominated, but the exact-mass matching step could never
# assign it, because the reference value was 2.14 Da away from anything the
# instrument could observe. The compound had to be confirmed by hand.
#
# This script prevents that class of failure from recurring. It recomputes the
# monoisotopic mass of every reference compound from its own molecular formula
# and reports any row where the tabulated value disagrees. It runs before any
# sample data is processed, so problems surface immediately rather than as
# unexplained absences in the final report.
# =============================================================================

validate_reference_library <- function(path, label,
                                       name_col, formula_col, mass_col,
                                       tol_mDa = 5) {

  if (!file.exists(path)) {
    warning("[", label, "] file not found: ", path)
    return(NULL)
  }

  lib <- read_csv(path, show_col_types = FALSE)

  # tolerate differences in header capitalisation / spacing between libraries
  pick <- function(target) {
    hit <- names(lib)[tolower(trimws(names(lib))) == tolower(target)]
    if (length(hit)) hit[1] else NA_character_
  }
  nc <- pick(name_col); fc <- pick(formula_col); mc <- pick(mass_col)

  if (is.na(fc) || is.na(mc)) {
    message(sprintf("  [%s] SKIPPED - no formula/mass column to cross-check", label))
    return(NULL)
  }

  chk <- lib %>%
    transmute(
      compound   = if (!is.na(nc)) .data[[nc]] else NA_character_,
      formula    = .data[[fc]],
      mass_given = suppressWarnings(as.numeric(.data[[mc]]))
    ) %>%
    filter(!is.na(formula), formula != "") %>%
    mutate(
      mass_calc = map_dbl(formula, formula_to_mono),
      delta_Da  = mass_given - mass_calc,
      delta_ppm = ifelse(is.na(mass_calc) | mass_calc == 0,
                         NA_real_, abs(delta_Da) / mass_calc * 1e6)
    )

  unparsed <- chk %>% filter(is.na(mass_calc))
  bad      <- chk %>% filter(!is.na(mass_calc), abs(delta_Da) > tol_mDa / 1000)

  message(sprintf("  [%s] %d entries | %d mass mismatches | %d unparseable formulae",
                  label, nrow(chk), nrow(bad), nrow(unparsed)))

  if (nrow(bad) > 0) {
    for (i in seq_len(min(nrow(bad), 10))) {
      message(sprintf("      %-42s %-12s given=%10.4f calc=%10.4f  d=%+8.4f Da",
                      substr(bad$compound[i], 1, 42), bad$formula[i],
                      bad$mass_given[i], bad$mass_calc[i], bad$delta_Da[i]))
    }
    if (nrow(bad) > 10) message(sprintf("      ... and %d more", nrow(bad) - 10))
  }

  bind_rows(
    bad      %>% mutate(issue = "mass does not match formula"),
    unparsed %>% mutate(issue = "formula could not be parsed")
  ) %>%
    mutate(library = label, .before = 1)
}


# -----------------------------------------------------------------------------
# Run the check across every reference table used by the pipeline.
# -----------------------------------------------------------------------------
run_library_validation <- function() {

  message("\n=== Reference library validation ===")

  report <- bind_rows(
    validate_reference_library(
      file.path(dir_data$wf1, "PAC_mass_spectral_library.csv"),
      "WF1 PAC spectral library", "Compound", "formula", "Exact mass"),

    validate_reference_library(
      file.path(dir_data$wf2, "possible_halo_pacs2.csv"),
      "WF2 halo-PAC exact masses", "possible_compounds", "Formula", "ExactMass"),

    validate_reference_library(
      file.path(dir_data$wf3, "NSO_PAC_exact_mass_reference.csv"),
      "WF3 NSO-PAC exact masses", "compound_name", "formula", "exact_mass"),

    validate_reference_library(
      file.path(dir_data$wf3, "HPAC_mass_spectral_library.csv"),
      "WF3 HPAC spectral library", "Compound", "formula", "Exact mass")
  )

  # Also check the KMD training standards' manually-assigned classes here, so
  # a mis-entered cell is caught at the start of the run rather than only when
  # the KMD stage happens to execute. run_kmd() repeats this same check
  # immediately before using the file, so it is still caught even if this
  # function is skipped.
  gt_file <- file.path(dir_data$wf2, "std_check_ground_truth_halo_class.csv")
  if (file.exists(gt_file)) {
    validate_standards_ground_truth(gt_file)
  }

  out <- file.path(dir_logs, paste0("library_validation_", Sys.Date(), ".csv"))

  if (is.null(report) || nrow(report) == 0) {
    message("  All reference libraries passed: every tabulated mass agrees with",
            " its formula to within 5 mDa.")
    write_csv(tibble(library = character(), compound = character(),
                     issue = character()), out)
  } else {
    write_csv(report, out)
    message("  ", nrow(report), " issue(s) written to ", basename(out))
    message("  ACTION: fix these entries before trusting affected assignments.")
  }

  invisible(report)
}


# -----------------------------------------------------------------------------
# validate_standards_ground_truth() — check the manually-assigned halogen
# class of each training standard against what its own molecular formula
# implies, using the same PAC_mass_spectral_library.csv already loaded by WF1.
# -----------------------------------------------------------------------------
# WHY THIS EXISTS
#
# The KMD stage (R/50_kmd_confidence.R) trains its regression models from a
# ground-truth file in which each standard's halogen_class has been entered
# by hand, deliberately bypassing WF2's own algorithmic guess (see that
# script's header comment for the reasoning). Manual entry across many rows
# is exactly the kind of step where a single mis-typed cell is easy to miss —
# this check caught "5,6-Br2-1,2-H2-Acy" entered as poly_or_mixed when its
# formula, C12H8Br2, is unambiguously di_Br, during preparation of this file.
#
# The check is deliberately independent of WF2: it derives the expected class
# purely by counting Cl and Br atoms in the compound's known formula, and
# compares that to whatever was typed into halogen_class. It has no opinion
# about isotope patterns or spectra at all.
# -----------------------------------------------------------------------------
expected_class_from_formula <- function(formula) {
  if (is.na(formula) || formula == "") return(NA_character_)
  m <- str_match_all(formula, "([A-Z][a-z]?)(\\d*)")[[1]]
  if (nrow(m) == 0) return(NA_character_)
  els <- m[, 2]
  cnt <- ifelse(m[, 3] == "", 1L, suppressWarnings(as.integer(m[, 3])))
  counts <- setNames(cnt, els)
  n_cl <- sum(counts[names(counts) == "Cl"])
  n_br <- sum(counts[names(counts) == "Br"])
  if (length(n_cl) == 0) n_cl <- 0
  if (length(n_br) == 0) n_br <- 0

  dplyr::case_when(
    n_cl == 1 && n_br == 0 ~ "mono_Cl",
    n_cl == 2 && n_br == 0 ~ "di_Cl",
    n_br == 1 && n_cl == 0 ~ "mono_Br",
    n_br == 2 && n_cl == 0 ~ "di_Br",
    (n_cl + n_br) > 0      ~ "poly_or_mixed",
    TRUE                   ~ NA_character_
  )
}

validate_standards_ground_truth <- function(
    ground_truth_file,
    library_file = file.path(dir_data$wf1, "PAC_mass_spectral_library.csv")) {

  message("\n=== Ground-truth halogen class check ===")

  if (!file.exists(library_file)) {
    warning("Ground-truth check: reference library not found, skipping: ",
            library_file)
    return(invisible(NULL))
  }

  lib <- read_csv(library_file, show_col_types = FALSE)
  lib_name_col <- names(lib)[tolower(trimws(names(lib))) == "acronyms"]
  lib_formula_col <- names(lib)[tolower(trimws(names(lib))) == "formula"]
  if (length(lib_name_col) == 0 || length(lib_formula_col) == 0) {
    warning("Ground-truth check: could not find Acronyms/formula columns in ",
            library_file, " - skipping")
    return(invisible(NULL))
  }

  expected <- lib %>%
    transmute(name = .data[[lib_name_col]],
             expected_class = map_chr(.data[[lib_formula_col]],
                                      expected_class_from_formula)) %>%
    filter(!is.na(expected_class))

  gt <- read_csv(ground_truth_file, show_col_types = FALSE) %>%
    transmute(name, entered_class = str_trim(str_extract(halogen_class, "^[^;]+")))

  mismatches <- gt %>%
    inner_join(expected, by = "name") %>%
    filter(entered_class != expected_class)

  not_in_library <- gt %>% anti_join(expected, by = "name")

  if (nrow(mismatches) == 0) {
    message("  All ", nrow(gt), " ground-truth entries agree with their ",
            "library formula.")
  } else {
    message("  ", nrow(mismatches), " entry(ies) DO NOT match their formula:")
    for (i in seq_len(nrow(mismatches))) {
      message(sprintf("    %-28s entered='%s'   formula implies='%s'",
                      mismatches$name[i], mismatches$entered_class[i],
                      mismatches$expected_class[i]))
    }
    message("  ACTION: fix these cells in ", basename(ground_truth_file),
            " before trusting the KMD models trained from it.")
  }

  if (nrow(not_in_library) > 0) {
    message("  Note: ", nrow(not_in_library), " entry(ies) not matched to a ",
            "library formula by name (renamed or not in this library) - not ",
            "checked: ", paste(not_in_library$name, collapse = ", "))
  }

  invisible(mismatches)
}

