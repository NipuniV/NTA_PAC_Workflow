# =============================================================================
# HIGH-CONFIDENCE HALO-PAC INTEGRATION — combine WF2 and KMD evidence
# -----------------------------------------------------------------------------
# PURPOSE
#   Produce a single ranked list, high_confidence_halo_PACs, by integrating the
#   two independent-ish lines of evidence the pipeline generates for a
#   halogenated peak:
#
#     WF2 (isotope pattern + exact-mass match)
#         - a valid Cl/Br isotope cluster is near-decisive evidence of halogen
#           content (fixed nuclide abundances; see WF2 header)
#         - an exact-mass match to the halo-PAC reference table proposes a
#           specific identity/class
#
#     KMD (homologous-series regression, theoretical + empirical training)
#         - tests whether the peak's Kendrick mass defect falls within the 95%
#           prediction interval of a halogen class's homologous-series line
#
# WHY KMD IS USED FOR CORROBORATION, NOT AS INDEPENDENT PROOF OF "HALO-PAC"
#
#   The KMD models are trained on the theoretical halo-PAC library
#   (possible_halo_pacs2.csv). Kendrick mass defect is a function of mass alone
#   and carries no ring-structure information, so "falls within the di-Br band"
#   really means "has a mass defect consistent with some di-brominated
#   C/H/Br composition" -- which a di-brominated NON-PAC could also satisfy.
#   KMD-in-band is therefore NOT clean independent evidence that a compound is a
#   halo-PAC; it substantially overlaps the "is this halogenated?" signal that
#   WF2's isotope step already establishes.
#
#   Treating KMD agreement as independent confirmation of halo-PAC identity
#   would be partly circular: it would reward a compound for resembling the very
#   set used to define the classes. This step therefore does NOT do that.
#   Instead KMD is used for the two things it CAN establish without circularity:
#
#     (1) CLASS CORROBORATION -- for a peak WF2 has already flagged as
#         halogenated and matched to a class, does the KMD evidence agree with
#         that same class? Agreement strengthens the class assignment.
#     (2) ANOMALY FLAGGING -- a halogenated peak that falls OUTSIDE all
#         prediction bands is behaving unlike any known halogen class: a flag
#         for manual review (possible novel structure, or a misassignment),
#         not a rejection.
#
# TIERED OUTPUT (each tier states what was independently established):
#
#   High          WF2 isotope-positive  AND  exact-mass match to reference
#                 AND  KMD in-band for the matched class
#                 -> three concordant lines; strongest assignment
#
#   Medium (a)    WF2 isotope-positive  AND  exact-mass match,
#                 but KMD out-of-band            [flag: class anomaly]
#                 -> identity proposed and halogen confirmed, but it does not
#                    sit on the expected homologous line for that class
#
#   Medium (b)    WF2 isotope-positive  AND  KMD in-band,
#                 but NO exact-mass match        [flag: possible novel halo-PAC]
#                 -> behaves like a known halogen class but matches nothing in
#                    the reference table: a candidate genuinely new compound
#
#   Low / review  WF2 isotope-positive only
#                 -> halogen content indicated, but neither corroborating line
#                    present; needs manual review
#
# INPUT   results_wf2/<sample>_all_halo_pacs.csv   (every halogen-positive peak)
#         results_wf2/<sample>_halo_pacs.csv       (subset with a reference match)
#         results_wf4/<sample>_kmd_ci_support.csv  (KMD per-compound result)
# OUTPUT  results_final/<sample>_high_confidence_halo_PACs.csv
# =============================================================================

run_high_confidence_halo_pacs <- function(sample_name) {

  message("\n=== HIGH-CONFIDENCE HALO-PAC INTEGRATION ===")

  p_all <- file.path(dir_res$wf2, paste0(sample_name, "_all_halo_pacs.csv"))
  p_lib <- file.path(dir_res$wf2, paste0(sample_name, "_halo_pacs.csv"))
  p_kmd <- file.path(dir_res$wf4, paste0(sample_name, "_kmd_ci_support.csv"))

  if (!file.exists(p_all)) {
    warning("HC halo-PAC: WF2 all-halo output not found - step skipped: ", p_all)
    return(invisible(NULL))
  }

  # ---- 6.1 WF2 evidence ----------------------------------------------------
  # all_halo   = every peak with a valid halogen isotope pattern (WF2 step)
  # matched    = the subset that ALSO matched the exact-mass reference table
  all_halo <- clean_names(read_csv(p_all, show_col_types = FALSE))

  matched <- if (file.exists(p_lib)) {
    clean_names(read_csv(p_lib, show_col_types = FALSE))
  } else {
    tibble(name = character(), rt1_s = numeric(), rt2_s = numeric())
  }

  # A peak is "isotope-positive" simply by being in all_halo. "reference-
  # matched" is being in the matched subset. Join keys are the peak coordinates
  # (name + both retention times), which uniquely identify a peak.
  keys <- c("name", "rt1_s", "rt2_s")

  wf2 <- all_halo %>%
    mutate(wf2_isotope_positive = TRUE) %>%
    left_join(
      matched %>%
        select(any_of(c(keys, "possible_halo_pac", "matched_formula",
                        "matched_via", "best_ppm_error"))) %>%
        mutate(wf2_reference_matched = TRUE),
      by = keys
    ) %>%
    mutate(wf2_reference_matched = coalesce(wf2_reference_matched, FALSE),
           # the halogen class WF2 assigned (first token, as elsewhere)
           wf2_class = str_trim(str_extract(halogen_class, "^[^;]+")))

  message(sprintf("  WF2 isotope-positive peaks : %d", nrow(wf2)))
  message(sprintf("  ...of which reference-matched: %d", sum(wf2$wf2_reference_matched)))

  # ---- 6.2 KMD evidence ----------------------------------------------------
  # From the KMD stage: for each peak, how many series placed it in-band, which
  # class the KMD evidence best supports, and the summary confidence label.
  if (file.exists(p_kmd)) {
    kmd <- clean_names(read_csv(p_kmd, show_col_types = FALSE)) %>%
      select(any_of(c(keys, "kmd_supported_class", "series_in_band",
                      "kmd_confidence", "extrapolation_flag"))) %>%
      distinct(across(all_of(keys)), .keep_all = TRUE) %>%
      mutate(kmd_in_band = coalesce(series_in_band, 0L) >= 1L)
  } else {
    warning("HC halo-PAC: KMD output not found - proceeding with WF2 only, ",
            "all KMD evidence treated as absent: ", p_kmd)
    kmd <- tibble(name = character(), rt1_s = numeric(), rt2_s = numeric(),
                  kmd_supported_class = character(), series_in_band = integer(),
                  kmd_confidence = character(), extrapolation_flag = character(),
                  kmd_in_band = logical())
  }

  # ---- 6.3 Integrate and assign tiers --------------------------------------
  integrated <- wf2 %>%
    left_join(kmd, by = keys) %>%
    mutate(
      kmd_in_band = coalesce(kmd_in_band, FALSE),
      # does the KMD-supported class agree with the class WF2 assigned?
      # only meaningful when both a WF2 class and an in-band KMD class exist
      class_agreement = case_when(
        !kmd_in_band                              ~ NA,
        is.na(wf2_class) | is.na(kmd_supported_class) ~ NA,
        wf2_class == kmd_supported_class          ~ TRUE,
        TRUE                                      ~ FALSE
      ),

      confidence_tier = case_when(
        wf2_reference_matched &  kmd_in_band ~ "High",
        wf2_reference_matched & !kmd_in_band ~ "Medium",
        !wf2_reference_matched &  kmd_in_band ~ "Medium",
        TRUE                                  ~ "Low / review"
      ),

      # a short, explicit statement of what each tier rests on + any flag
      evidence_basis = case_when(
        wf2_reference_matched &  kmd_in_band ~
          "isotope + reference match + KMD in-band",
        wf2_reference_matched & !kmd_in_band ~
          "isotope + reference match; KMD out-of-band [class anomaly - review]",
        !wf2_reference_matched &  kmd_in_band ~
          "isotope + KMD in-band; no reference match [possible novel halo-PAC]",
        TRUE ~
          "isotope pattern only [no corroboration - review]"
      ),

      # promote the class-disagreement case to an explicit flag even within High
      review_flag = case_when(
        confidence_tier == "High" & class_agreement == FALSE ~
          "KMD supports a DIFFERENT class than WF2 - check class assignment",
        !is.na(extrapolation_flag) & extrapolation_flag != "" ~ extrapolation_flag,
        TRUE ~ NA_character_
      )
    )

  # ---- 6.4 Assemble the output table ---------------------------------------
  tier_rank <- c("High" = 1, "Medium" = 2, "Low / review" = 3)

  final <- integrated %>%
    transmute(
      name, rt1_s, rt2_s,
      confidence_tier,
      evidence_basis,
      review_flag,
      wf2_class,
      wf2_isotope_positive,
      wf2_reference_matched,
      possible_halo_pac,
      matched_formula,
      best_ppm_error,
      kmd_supported_class,
      kmd_in_band,
      series_in_band,
      class_agreement,
      kmd_confidence,
      aromaticity_index,
      detected_isotope_peaks
    ) %>%
    arrange(tier_rank[confidence_tier], name, rt1_s, rt2_s)

  out <- file.path(dir_final, paste0(sample_name, "_high_confidence_halo_PACs.csv"))
  write_csv(final, out)

  message(sprintf("  Integrated peaks : %d", nrow(final)))
  print(count(final, confidence_tier, name = "n") %>%
          arrange(match(confidence_tier, names(tier_rank))))
  n_flagged <- sum(!is.na(final$review_flag))
  if (n_flagged > 0) {
    message(sprintf("  Peaks with a review flag: %d", n_flagged))
  }
  message("  Saved: results_final/", basename(out))

  invisible(final)
}
