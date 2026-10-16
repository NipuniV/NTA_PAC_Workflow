# =============================================================================
# WORKFLOW 3 — Screening for N-, S- and O-heterocyclic PACs (HPACs)
# -----------------------------------------------------------------------------
# PURPOSE
#   Identify polycyclic aromatic compounds in which one or more ring carbons is
#   replaced by nitrogen, sulfur or oxygen.
#
# WHY THIS IS HARDER THAN WF2 — and why the workflow is built differently
#
#   WF2 can rely on a physical constant: chlorine and bromine isotope ratios are
#   fixed properties of the nuclide, so one clean isotope cluster is close to
#   proof of halogen content.
#
#   N, S and O offer no equivalent signature under electron ionisation. Their
#   contribution to the mass defect is small relative to carbon and hydrogen,
#   and the neutral losses they produce (HCN, CO, H2S) are also produced by many
#   ordinary non-aromatic organics. Heteroatom screening is therefore inherently
#   INFERENTIAL: several individually weak indicators combined, rather than one
#   decisive one.
#
#   WF3 is consequently designed for HIGH RECALL, not high precision. An HPAC
#   discarded here can never be recovered later, whereas a false candidate can
#   be removed downstream by WF4. Its output is a candidate list, never a set of
#   identifications.
#
# INPUT   data_raw/<sample>.csv                  (shared with WF1 and WF2)
#         data_wf3/NSO_PAC_exact_mass_reference.csv
#         data_wf3/HPAC_mass_spectral_library.csv
# OUTPUT  results_wf3/NSO_PAC_<sample>_final.csv   (all screened features)
#         results_wf3/<sample>_hpacs.csv           (high-confidence subset)
#
# NOTE: this workflow deliberately STOPS after writing the CSVs. The KMD and
#       contour plots that previously followed are not produced here; the
#       consolidated KMD analysis lives in R/50_kmd_confidence.R.
# =============================================================================


# -----------------------------------------------------------------------------
# check_all_ion_losses() — test every ordered pair of ions in a spectrum for the
# diagnostic neutral losses.
#
# The ppm error is referenced to the LOWER (fragment) ion because that is the
# ion actually measured by the instrument, and absolute mass uncertainty scales
# with measured mass. This gives a window that tightens at low mass and relaxes
# at high mass, in proportion to real instrument performance — unlike a fixed
# Da window, which is too permissive at low mass and too strict at high mass.
# -----------------------------------------------------------------------------
check_all_ion_losses <- function(mass_list, tol_ppm = wf3_frag_tol_ppm) {
  mass_v <- sort(unique(mass_list))
  n <- length(mass_v)
  as.list(vapply(NEUTRAL_LOSSES, function(loss) {
    found <- FALSE
    for (i in seq_len(n)) {
      for (j in seq_len(i - 1L)) {
        if (abs((mass_v[i] - mass_v[j]) - loss) / mass_v[j] * 1e6 <= tol_ppm) {
          found <- TRUE; break
        }
      }
      if (found) break
    }
    found
  }, logical(1)))
}


# -----------------------------------------------------------------------------
# match_any_mass() — match ANY fragment mass against the NSO reference table.
#
# IMPORTANT — read before interpreting the ref_compound column:
# This deliberately scans every fragment, not only the base peak, because in EI
# spectra the base peak is frequently a fragment rather than the molecular ion.
# The cost is a large number of comparisons per feature (typically ~50 fragments
# x 154 reference entries), which makes coincidental hits possible. The
# reference table also contains 56 alkyl homologues spaced by CH2, and EI
# fragment series are spaced the same way, which further inflates the chance of
# an accidental match.
#
# A reference match is therefore CANDIDATE GENERATION, not evidence of identity.
# -----------------------------------------------------------------------------
match_any_mass <- function(mass_list, ref_df, tol_ppm) {
  best_ppm <- Inf; best_idx <- NA_integer_
  for (obs_mass in mass_list) {
    ppm <- abs(obs_mass - ref_df$exact_mass) / obs_mass * 1e6
    i   <- which.min(ppm)
    if (length(i) && ppm[i] < best_ppm) { best_ppm <- ppm[i]; best_idx <- i }
  }
  if (is.na(best_idx) || best_ppm > tol_ppm) {
    return(list(ref_match = FALSE, ref_compound = NA_character_,
                ref_formula = NA_character_, ref_type = NA_character_,
                ref_ppm = NA_real_))
  }
  list(ref_match    = TRUE,
       ref_compound = ref_df$compound_name[best_idx],
       ref_formula  = ref_df$formula[best_idx],
       ref_type     = ref_df$compound_type[best_idx],
       ref_ppm      = best_ppm)
}


run_wf3 <- function(sample_name) {

  message("\n=== WORKFLOW 3 : NSO heterocyclic PAC screen ===")
  message("  Scoring scheme: ", wf3_scoring)

  # ---- 3.1 Load ------------------------------------------------------------
  sample_raw <- read_spectra(raw_path(sample_name)) %>%
    mutate(across(2:ncol(.), as.numeric))
  nso_ref <- read_csv(file.path(dir_data$wf3, "NSO_PAC_exact_mass_reference.csv"),
                      show_col_types = FALSE) %>%
    mutate(exact_mass = as.numeric(exact_mass))
  library_raw <- read_spectra(file.path(dir_data$wf3, "HPAC_mass_spectral_library.csv"))

  # ---- 3.2 Noise and column-bleed removal ---------------------------------
  filtered_base <- sample_raw %>%
    filter(x1st_dimension_time_s >= wf3_min_rt1,
           base_mass             >= wf3_min_base_mass)
  message(sprintf("  After solvent/low-mass cut : %d", n_distinct(filtered_base$name)))

  sample_clean <- remove_column_bleed(filtered_base, wf3_bleed_ppm)
  message(sprintf("  After bleed removal        : %d", n_distinct(sample_clean$name)))

  # ---- 3.3 Reshape: ALL ions retained -------------------------------------
  # Unlike WF1, no top-N truncation: the neutral-loss test needs every ion pair.
  samp_long <- reshape_ms_long(sample_clean, "name", top_n = NULL) %>%
    rename(frag_mass = m, frag_int = ri)

  samp_features <- samp_long %>%
    group_by(name, base_mass, x1st_dimension_time_s, x2nd_dimension_time_s) %>%
    summarise(all_masses  = list(sort(unique(frag_mass))),
              n_fragments = n_distinct(frag_mass),
              .groups = "drop")

  # ---- 3.4 Reference exact-mass matching ----------------------------------
  samp_features <- samp_features %>%
    mutate(ref_result   = map(all_masses, match_any_mass,
                              ref_df = nso_ref, tol_ppm = wf3_mass_tol_ppm),
           ref_match    = map_lgl(ref_result, "ref_match"),
           ref_compound = map_chr(ref_result, "ref_compound"),
           ref_formula  = map_chr(ref_result, "ref_formula"),
           ref_type     = map_chr(ref_result, "ref_type"),
           ref_ppm      = map_dbl(ref_result, "ref_ppm")) %>%
    select(-ref_result)
  message(sprintf("  Features with a reference match : %d of %d",
                  sum(samp_features$ref_match), nrow(samp_features)))

  # ---- 3.5 Kendrick mass defect (annotation only, not a filter) -----------
  samp_kmd <- samp_features %>%
    mutate(kendrick_mass = base_mass * (14 / 14.01565),
           kmd           = round(kendrick_mass) - kendrick_mass,
           kmd_flag      = abs(kmd) > 0.05)

  # ---- 3.6 Isotope checks --------------------------------------------------
  check_isotopes <- function(mass_list, parent_mass, tol = wf3_iso_tol) {
    list(has_M1 = any(abs(mass_list - (parent_mass + 1.00335)) <= tol),
         has_M2 = any(abs(mass_list - (parent_mass + 2.00671)) <= tol))
  }
  samp_iso <- samp_kmd %>%
    mutate(iso_check = map2(all_masses, base_mass, check_isotopes),
           has_M1    = map_lgl(iso_check, "has_M1"),
           has_M2    = map_lgl(iso_check, "has_M2")) %>%
    select(-iso_check)

  # ---- 3.7 Aromaticity (degree of unsaturation proxy) ---------------------
  samp_dbe <- samp_iso %>%
    mutate(odd_nominal  = (round(base_mass) %% 2) == 1,
           dbe_aromatic = (base_mass / 26) >= 4)

  # ---- 3.8 Diagnostic neutral losses --------------------------------------
  samp_frag <- samp_dbe %>%
    mutate(frag_check  = map(all_masses, check_all_ion_losses),
           frag_N_loss = map_lgl(frag_check, "HCN"),
           frag_S_loss = map_lgl(frag_check, "S"),
           frag_O_loss = map_lgl(frag_check, "CO"),
           frag_H2S    = map_lgl(frag_check, "H2S"),
           frag_OH     = map_lgl(frag_check, "OH"),
           frag_evidence = case_when(
             frag_N_loss & frag_S_loss & frag_O_loss ~ "N+S+O",
             frag_N_loss & frag_S_loss               ~ "N+S",
             frag_N_loss & frag_O_loss               ~ "N+O",
             frag_S_loss & frag_O_loss               ~ "S+O",
             frag_N_loss                             ~ "N",
             frag_S_loss                             ~ "S",
             frag_O_loss                             ~ "O",
             TRUE                                    ~ "None"),
           sulfur_iso_flag = has_M2 & (frag_S_loss | frag_H2S)) %>%
    select(-frag_check)

  # ---- 3.9 Evidence scoring ------------------------------------------------
  #
  # Three schemes are available; set wf3_scoring in R/00_setup.R.
  #
  # "original" — the scheme used for the thesis and submitted manuscript.
  #     max 8 points:
  #       +2 kmd_flag, +1 has_M1, +1 has_M2, +1 sulfur_iso_flag,
  #       +1 any neutral loss, +1 frag_N_loss, +1 dbe_aromatic
  #
  # "recetox" — re-weighted using the RECETOX Exposome benchmark (386 spectra,
  #     23 true HPACs). Each criterion was measured for DISCRIMINATION, defined
  #     as P(fires | true HPAC) - P(fires | not HPAC):
  #
  #       criterion            fires overall   discrimination
  #       kmd_flag (2 pts)         93.9%           +0.018
  #       dbe_aromatic             99.2%           +0.009
  #       has_M1                   93.4%           +0.071
  #       has_M2                   16.0%           -0.031
  #       sulfur_iso_flag           1.9%           -0.021
  #       any neutral loss         79.9%           +0.122
  #       frag_N_loss (HCN)        35.0%           +0.137
  #
  #     Five of the eight original points therefore came from criteria that fire
  #     for essentially every feature and carry almost no class information; two
  #     were slightly negative. A background compound with no heteroatom evidence
  #     scored 5/8 by default, which is exactly why 293 of 386 RECETOX features
  #     landed in the top tier.
  #
  # "minimal" (DEFAULT) — the manuscript-stage scheme. Keeps ONLY the two
  #     criteria with real, defensible chemical discriminatory power and drops
  #     the four weak ones entirely (KMD, M+1, M+2, aromatic DBE), rather than
  #     down-weighting them. The rationale for each retained/dropped criterion:
  #
  #       KEPT:
  #       +1 ref_match (fragment-ion match to PAC/HPAC reference DB, +/-5 ppm)
  #            Under EI, PACs fragment extensively before the molecular ion is
  #            recorded, so comparing ALL detected fragment masses (not just the
  #            base peak) against known PAC/HPAC reference masses still catches
  #            diagnostic substructure fragments. Anchors a candidate to known
  #            PAC-like fragmentation. Discriminatory power: moderate.
  #       +2 any diagnostic neutral loss (HCN / CO / S / H2S)
  #            Heteroatom rings fragment along heteroatom-specific pathways:
  #            N-rings expel HCN, O-rings expel CO, S-rings expel S/H2S. The
  #            identity of the lost neutral maps directly onto which heteroatom
  #            is present. The ONLY criterion shown to meaningfully separate
  #            HPACs from background. Discriminatory power: strong.
  #
  #       DROPPED (weak / near-universal / negative discrimination):
  #         KMD (CH2 base)   fired on 93.9% of ALL features - flags "alkylated"
  #                          more than "heteroatom-containing"
  #         M+1 (13C)        present in essentially every carbon-containing
  #                          compound; very weak
  #         M+2 (34S)        genuine sulfur pointer but shared with any
  #                          S-compound; slightly negative discrimination
  #         aromatic DBE     fired on 99.2% of features; true of homocyclic
  #                          PAHs too, so does not isolate HPACs
  #
  #     max 3 points:  +1 ref_match, +2 any diagnostic neutral loss
  #
  # Both "recetox" and "minimal" change WF3 output relative to the published
  # version. Use "original" to reproduce thesis numbers.

  if (identical(wf3_scoring, "original")) {
    samp_scored <- samp_frag %>%
      mutate(
        evidence_score =
          (kmd_flag * 2L) + (has_M1 * 1L) + (has_M2 * 1L) +
          (sulfur_iso_flag * 1L) + ((frag_evidence != "None") * 1L) +
          (frag_N_loss * 1L) + (dbe_aromatic * 1L),
        max_score = 8L
      )
  } else if (identical(wf3_scoring, "recetox")) {
    samp_scored <- samp_frag %>%
      # mandatory gate: a fused aromatic system is a prerequisite, so features
      # failing it score 0 rather than accumulating points from weak evidence
      mutate(
        evidence_score = if_else(
          dbe_aromatic,
          (frag_N_loss * 3L) +
            ((frag_S_loss | frag_H2S) * 2L) +
            (frag_O_loss * 2L) +
            (sulfur_iso_flag * 1L) +
            (has_M1 * 1L),
          0L),
        max_score = 9L
      )
  } else {
    # "minimal" — only the two criteria with genuine discriminatory power.
    # ref_match is the fragment-ion match to the PAC/HPAC reference DB (already
    # computed in step 3.4). "any diagnostic neutral loss" is TRUE whenever an
    # HCN, CO, S or H2S loss was detected (frag_evidence != "None").
    samp_scored <- samp_frag %>%
      mutate(
        evidence_score = (ref_match * 1L) + ((frag_evidence != "None") * 2L),
        max_score = 3L
      )
  }

  samp_scored <- samp_scored %>%
    mutate(
      score_fraction = evidence_score / max_score,
      screening_tier = case_when(
        score_fraction >= wf3_tier_high ~ "High - likely NSO-PAC",
        score_fraction >= wf3_tier_med  ~ "Medium - probable NSO-PAC",
        evidence_score >  0             ~ "Low - weak evidence",
        TRUE                            ~ "None - likely pure PAH / noise"
      )
    )

  # ---- 3.10 Optional annotation against the in-house HPAC library ----------
  # Provides an orthogonal line of evidence (real measured spectrum + real
  # retention time) rather than a purely computational inference.
  parse_formula <- function(formula) {
    get_count <- function(el) {
      m <- regmatches(formula, regexpr(paste0(el, "([0-9]*)"), formula))
      if (length(m) == 0 || m == "") return(0L)
      n <- sub(el, "", m); if (n == "") 1L else as.integer(n)
    }
    list(C = get_count("C"), H = get_count("H"),
         N = get_count("N"), O = get_count("O"), S = get_count("S"))
  }
  compute_dbe <- function(formula) {
    f <- parse_formula(formula); (2 * f$C + 2 + f$N - f$H) / 2
  }

  library_nso <- library_raw %>%
    mutate(dbe   = map_dbl(formula, possibly(compute_dbe, NA_real_)),
           has_N = str_detect(formula, "N"),
           has_S = str_detect(formula, "S"),
           has_O = str_detect(formula, "O"),
           is_nso = has_N | has_S | has_O) %>%
    filter(!is.na(dbe), dbe >= 4, is_nso)
  message(sprintf("  HPAC library entries usable : %d", nrow(library_nso)))

  lib_top4  <- reshape_ms_long(library_nso, "acronyms", top_n = 4)
  samp_top4 <- samp_long %>% group_by(name) %>%
    slice_max(rel_int_pct, n = 4, with_ties = FALSE) %>% ungroup()

  lib_summary <- library_nso %>%
    select(acronyms, formula, dbe, has_N, has_S, has_O,
           base_mass, x1st_dimension_time_s, x2nd_dimension_time_s) %>%
    rename(base_mass_lib = base_mass, t1_lib = x1st_dimension_time_s,
           t2_lib = x2nd_dimension_time_s)

  samp_summary <- samp_scored %>%
    select(name, base_mass, x1st_dimension_time_s, x2nd_dimension_time_s) %>%
    rename(t1_samp = x1st_dimension_time_s, t2_samp = x2nd_dimension_time_s)

  ion_overlap <- samp_top4 %>%
    select(name, frag_mass) %>%
    cross_join(lib_top4 %>% select(acronyms, m)) %>%
    filter(abs(frag_mass - m) < 0.005) %>%
    group_by(name, acronyms) %>%
    summarise(ions_matched_count = n(), .groups = "drop")

  lib_hits <- samp_summary %>%
    cross_join(lib_summary) %>%
    mutate(ppm_diff = abs(base_mass - base_mass_lib) / base_mass_lib * 1e6) %>%
    filter(ppm_diff <= wf3_mass_tol_ppm) %>%
    left_join(ion_overlap, by = c("name", "acronyms")) %>%
    replace_na(list(ions_matched_count = 0L)) %>%
    filter(ions_matched_count >= 2) %>%
    mutate(t1_dist_s = abs(t1_samp - t1_lib),
           t2_dist_s = abs(t2_samp - t2_lib),
           rt_match  = (t1_dist_s <= t1_tol & t2_dist_s <= t2_tol),
           lib_confidence = case_when(
             rt_match & ions_matched_count >= 3 ~ "Library L1 - RT + 3 ions",
             rt_match & ions_matched_count >= 2 ~ "Library L2 - RT + 2 ions",
             ions_matched_count >= 3            ~ "Library L3 - mass/ions (isomer?)",
             TRUE                               ~ "Library L4 - mass match only"))

  best_hits <- lib_hits %>%
    arrange(name, desc(ions_matched_count), t1_dist_s) %>%
    group_by(name) %>% slice(1) %>% ungroup() %>%
    select(name, acronyms, formula, dbe, has_N, has_S, has_O,
           ppm_diff, ions_matched_count, rt_match, t1_dist_s, t2_dist_s, lib_confidence)

  all_lib_hits <- lib_hits %>%
    arrange(name, desc(ions_matched_count), t1_dist_s) %>%
    group_by(name) %>%
    summarise(all_lib_candidates = paste(acronyms, collapse = " | "), .groups = "drop")

  # ---- 3.11 Final table and export ----------------------------------------
  final_table <- samp_scored %>%
    select(-all_masses) %>%
    left_join(best_hits,    by = "name") %>%
    left_join(all_lib_hits, by = "name") %>%
    arrange(desc(evidence_score), screening_tier) %>%
    select(name, x1st_dimension_time_s, x2nd_dimension_time_s, base_mass,
           ref_compound, ref_formula, ref_type, ref_ppm,
           evidence_score, max_score, score_fraction, screening_tier,
           kmd, kmd_flag, has_M1, has_M2, sulfur_iso_flag, frag_evidence,
           frag_N_loss, frag_S_loss, frag_O_loss, frag_H2S, frag_OH,
           dbe_aromatic, odd_nominal, n_fragments,
           lib_confidence, acronyms, formula, dbe, has_N, has_S, has_O,
           ppm_diff, ions_matched_count, rt_match, t1_dist_s, t2_dist_s,
           all_lib_candidates)

  out_full <- file.path(dir_res$wf3, paste0("NSO_PAC_", sample_name, "_final.csv"))
  write_csv(final_table, out_full)
  message(sprintf("  Features reported : %d", nrow(final_table)))
  print(count(final_table, screening_tier, name = "n"))
  message("  Saved: ", basename(out_full))

  # ---- 3.12 High-confidence subset ----------------------------------------
  high_conf <- final_table %>% filter(score_fraction >= wf3_tier_high)
  out_high  <- file.path(dir_res$wf3, paste0(sample_name, "_hpacs.csv"))
  write_csv(high_conf, out_high)
  message(sprintf("  High-confidence   : %d", nrow(high_conf)))
  message("  Saved: ", basename(out_high))

  stage_file(out_high, dir_data$wf4)

  # WF3 ends here by design. KMD plotting is handled in R/50_kmd_confidence.R.
  invisible(list(full = final_table, high = high_conf))
}
