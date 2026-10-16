# =============================================================================
# run_all.R — Single entry point for the complete PAC non-targeted workflow
# -----------------------------------------------------------------------------
# Runs, in order:
#
#   0. Reference-library validation   (checks every tabulated mass against its
#                                      own formula BEFORE any data is processed)
#   1. WF1  library matching and confidence levels
#   2. WF2  halogenated PAC isotope filter
#   3. WF3  NSO heterocyclic PAC screen
#   4. WF4  integration of all three evidence streams
#   5. KMD  homologous-series confidence + the three final plots
#
# WF1, WF2 and WF3 are INDEPENDENT of one another and could be run in any order;
# only WF4 and the KMD stage depend on earlier outputs. Each is kept in its own
# script so that individual workflows can still be run, cited and described
# separately.
#
# HOW TO USE
#   1. Put your ChromaTOF export in data_wf1/, data_wf2/ and data_wf3/
#   2. Set `sample_name` in R/00_setup.R
#   3. source("run_all.R")
#
# To run one workflow on its own:
#   source("R/00_setup.R"); source("R/01_helpers.R"); source("R/20_wf2_halogen_filter.R")
#   run_wf2(sample_name)
# =============================================================================

t_start <- Sys.time()

source(here::here("R", "00_setup.R"))
source(here::here("R", "01_helpers.R"))
source(here::here("R", "02_validate_libraries.R"))
source(here::here("R", "10_wf1_library_match.R"))
source(here::here("R", "20_wf2_halogen_filter.R"))
source(here::here("R", "30_wf3_hpac_filter.R"))
source(here::here("R", "40_wf4_integration.R"))
source(here::here("R", "50_kmd_confidence.R"))
source(here::here("R", "60_high_confidence_halo_pacs.R"))

message("\n", strrep("=", 70))
message("PAC NON-TARGETED ANALYSIS PIPELINE")
message("Sample: ", sample_name)
message(strrep("=", 70))

# ---- Step 0: validate reference libraries before touching any sample data ---
# A wrong tabulated mass causes a compound to disappear silently rather than
# raise an error, so this check runs first and reports problems up front.
run_library_validation()

# ---- Steps 1-3: the three independent filters -------------------------------
# Wrapped so that a failure in one workflow does not prevent the others from
# running; WF4 handles missing inputs gracefully.
safely_run <- function(f, label, ...) {
  out <- tryCatch(f(...), error = function(e) {
    msg <- conditionMessage(e)
    # Print prominently AND record to disk. A stage that fails must never be
    # mistaken for a stage that simply found nothing.
    message("\n", strrep("!", 70))
    message("  STAGE FAILED: ", label)
    message("  ", msg)
    message(strrep("!", 70), "\n")
    cat(sprintf("[%s] %s FAILED: %s\n", format(Sys.time()), label, msg),
        file = file.path(dir_logs, paste0("run_errors_", Sys.Date(), ".log")),
        append = TRUE)
    NULL
  })
  invisible(out)
}

safely_run(run_wf1, "WF1", sample_name)
safely_run(run_wf2, "WF2", sample_name)
safely_run(run_wf3, "WF3", sample_name)

# ---- Step 4: integrate -------------------------------------------------------
safely_run(run_wf4, "WF4", sample_name)

# ---- Step 5: KMD model fitting + standards-only diagnostic plots ------------
# Independent of any sample. Fits the 12 regression models and saves
# kmd_final/standards_kmd_*.png (fitted lines + training data, no sample
# overlaid) BEFORE any sample-specific data is processed against them.
kmd_models <- safely_run(build_kmd_models, "KMD models")

# ---- Step 6: test the sample against those models, save sample plots --------
safely_run(run_kmd, "KMD", sample_name, models = kmd_models)

# ---- Step 7: integrate WF2 + KMD into the tiered high-confidence list -------
# Combines the two independent-ish halogen evidence streams into
# results_final/<sample>_high_confidence_halo_PACs.csv. KMD is used for class
# corroboration and anomaly flagging, not as independent proof of halo-PAC
# identity (see the header of R/60_high_confidence_halo_pacs.R for why).
safely_run(run_high_confidence_halo_pacs, "HC halo-PACs", sample_name)

# ---- Summary ----------------------------------------------------------------
elapsed <- round(as.numeric(difftime(Sys.time(), t_start, units = "secs")), 1)

message("\n", strrep("=", 70))
message("RUN COMPLETE  (", elapsed, " s)")
message(strrep("=", 70))
message("Primary outputs:")
message("  results_final/", sample_name, "_final_PAC_report.csv")
message("  results_final/", sample_name, "_high_confidence_halo_PACs.csv")
message("  kmd_final/standards_kmd_HCl.png  (model-only diagnostic, no sample data)")
message("  kmd_final/standards_kmd_Br2.png  (model-only diagnostic, no sample data)")
message("  kmd_final/standards_kmd_CH2.png  (model-only diagnostic, no sample data)")
message("  kmd_final/",     sample_name, "_kmd_HCl.png")
message("  kmd_final/",     sample_name, "_kmd_Br2.png")
message("  kmd_final/",     sample_name, "_kmd_CH2.png")
message("")
message("Intermediate outputs retained in results_wf1..4/ for traceability.")
message("Library validation report in logs/")
