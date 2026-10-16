# =============================================================================
# 00_setup.R — Global configuration for the integrated PAC NTA workflow
# -----------------------------------------------------------------------------
# This file is sourced first by run_all.R and by every individual workflow
# script. It defines:
#   (a) the packages used throughout,
#   (b) the project folder layout,
#   (c) every numerical tolerance / threshold used by WF1-WF4.
#
# ALL tunable parameters live here and NOWHERE else. If you need to change a
# tolerance, change it here once and every workflow picks it up. This is the
# single most important difference from the original separate scripts, where
# the same threshold appeared in several files and could drift out of sync.
# =============================================================================

# -----------------------------------------------------------------------------
# PACKAGES
# -----------------------------------------------------------------------------
# Any package that is not already installed is installed automatically, so the
# project runs on a clean machine without the user having to prepare anything
# first. This is what makes the pipeline genuinely transferable to another
# laboratory or reviewer.
required_packages <- c(
  "tidyverse",   # dplyr / tidyr / purrr / stringr / ggplot2
  "here",        # project-relative file paths (no absolute paths anywhere)
  "janitor",     # clean_names() for consistent snake_case columns
  "ggrepel",     # non-overlapping point labels on the KMD plots
  "knitr",       # R Markdown rendering
  "rmarkdown"    # R Markdown rendering
)

missing_packages <- required_packages[
  !(required_packages %in% rownames(installed.packages()))
]

if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(here)
  library(janitor)
  library(ggrepel)
})

# -----------------------------------------------------------------------------
# USER INPUT — the only line you normally need to edit
# -----------------------------------------------------------------------------
# Name of the ChromaTOF export to process, WITHOUT the .csv extension.
# The same base name is used for every intermediate and final output file so
# that results from different samples never overwrite one another.
#
# The raw file lives in ONE place only: data_raw/. WF1, WF2 and WF3 all read it
# from there. They are parallel filters, not a linear chain — only WF4 depends
# on the outputs of the other three.
sample_name <- "blank"

# -----------------------------------------------------------------------------
# FOLDER LAYOUT
# -----------------------------------------------------------------------------
# data_raw     : the ChromaTOF export, stored once and read by WF1, WF2 and WF3
# data_wfN     : the reference tables that workflow N needs (no sample data)
# results_wfN  : every file workflow N writes out
# results_final: the single integrated report from WF4
# kmd_final    : the three combined KMD plots (H/Cl, Br2, CH2)
# logs         : reference-library validation reports and run logs
dir_raw <- here("data_raw")   # single home for the ChromaTOF export

dir_data <- list(
  wf1 = here("data_wf1"),
  wf2 = here("data_wf2"),
  wf3 = here("data_wf3"),
  wf4 = here("data_wf4")
)
dir_res <- list(
  wf1 = here("results_wf1"),
  wf2 = here("results_wf2"),
  wf3 = here("results_wf3"),
  wf4 = here("results_wf4")
)
dir_final <- here("results_final")
dir_kmd   <- here("kmd_final")
dir_logs  <- here("logs")

# create any folder that does not yet exist (safe to re-run)
for (p in c(dir_raw, unlist(dir_data), unlist(dir_res),
            dir_final, dir_kmd, dir_logs)) {
  dir.create(p, showWarnings = FALSE, recursive = TRUE)
}

# raw_path() — full path to the ChromaTOF export for the current sample.
# Defined here so that all three workflows resolve the same single file.
raw_path <- function(sample_name) file.path(dir_raw, paste0(sample_name, ".csv"))

# -----------------------------------------------------------------------------
# SHARED TOLERANCES
# -----------------------------------------------------------------------------
# Retention-time windows used to decide whether a sample peak sits in the same
# 2D-GC position as a library standard. t1 = first dimension (volatility),
# t2 = second dimension (polarity/shape).
t1_tol <- 5      # seconds, 1st dimension
t2_tol <- 0.5    # seconds, 2nd dimension

# Isomer window: peaks sharing a base mass but eluting further apart than this
# are treated as unrelated rather than as isomers of the same compound.
t1_isomer_window <- 180  # seconds (3 minutes)

# Mass-accuracy tolerance for matching an observed base mass to a library
# standard. 5 ppm is appropriate for the GCxGC-HRTOF-MS platform (R = 25,000).
mass_tol_ppm <- 5

# Minimum number of the top-4 ions that must match for a library assignment
# to be considered at all (WF1).
min_ions_matched <- 1

# -----------------------------------------------------------------------------
# WF2 (halogen filter) parameters
# -----------------------------------------------------------------------------
wf2_min_rt1        <- 850     # s   drop solvent front / early-eluting light matrix
wf2_min_base_mass  <- 100     # Da  drop low-mass chemical noise
wf2_tol_mass       <- 0.002   # Da  window for pairing M with M+2
wf2_min_b_pct_of_a <- 0.40    # minimum M/(M+2) intensity ratio for a valid pair
wf2_iso_offset     <- 1.9970  # Da  nominal M -> M+2 spacing (Cl and Br alike)
wf2_ppm_match      <- 20      # ppm tolerance when matching to the halo-PAC table
wf2_bleed_ppm      <- 10      # ppm tolerance for column-bleed removal

# Aromaticity index acceptance window. The fold threshold is 0.8 (not 0.5) and
# the window is TWO-SIDED: a candidate is kept only when
# wf2_ai_min < AI < wf2_ai_max.
#
# Bounds re-optimized against the 824-compound AI validation set, with
# halo-PACs as the sole target class (WF2's actual objective). Findings:
#   - Lower bound tightened from -0.05 to -0.008: nothing legitimate sits below
#     -0.008 (the most negative real halo-PAC is -0.008), so the old -0.05 floor
#     was dead range that only admitted noise.
#   - Upper bound kept at 0.005: this gave the best F1 (0.703) of any window for
#     the halo-PAC target. Widening it does NOT help - 38 of 300 halo-PACs are
#     heavily brominated (avg 3.5 Br) and sit at AI 0.009-0.030, overlapping the
#     aliphatic interferent region, because Br's large negative mass defect
#     pushes AI upward. No window separates those poly-Br species from alkanes,
#     so widening the ceiling to catch them collapses specificity (0.65 -> 0.02)
#     while admitting every alkane. Those 38 are instead recovered downstream by
#     the isotope-pattern step, which sees their very strong, unambiguous Br
#     clusters - the AI gate is a coarse aromaticity pre-filter, not the halogen
#     detector.
wf2_ai_fold <- 0.8
wf2_ai_max  <- 0.005    # upper bound: best F1 for halo-PACs; see note above
wf2_ai_min  <- -0.008   # lower bound: tightened from -0.05 (dead range removed)

# Polysiloxane column-bleed signature masses (GC septa/column bleed).
bleed_masses <- c(73.0474, 207.0325, 221.0482, 253.0200,
                  267.0536, 281.0693, 355.0851, 429.1008)

# -----------------------------------------------------------------------------
# WF3 (NSO heterocyclic PAC filter) parameters
# -----------------------------------------------------------------------------
wf3_min_rt1       <- 800    # s
wf3_min_base_mass <- 100    # Da
wf3_mass_tol_ppm  <- 5      # ppm, reference exact-mass matching
wf3_iso_tol       <- 0.005  # Da, M+1 / M+2 isotope detection
wf3_frag_tol_ppm  <- 10     # ppm, diagnostic neutral-loss matching
wf3_bleed_ppm     <- 5      # ppm, column-bleed removal

# --- Evidence scoring scheme -------------------------------------------------
# "original"  reproduces the scoring used for the thesis / submitted manuscript.
# "recetox"   applies the re-weighting derived from the RECETOX benchmark, in
#             which four of the seven original criteria were shown to fire for
#             >90% of ALL features and therefore carried almost no information
#             (see R/30_wf3_hpac_filter.R for the measured discrimination of
#             each criterion). The re-weighted scheme moves the available points
#             onto the diagnostic neutral losses, which were the only criteria
#             that separated true HPACs from background.
#
# IMPORTANT: neither "recetox" nor "minimal" reproduces the numbers reported in
# the thesis. Use "original" whenever you need to regenerate published results.
#
# "minimal" (DEFAULT, manuscript stage) keeps ONLY the two criteria with
# genuine discriminatory power - fragment-ion match to the PAC/HPAC reference
# DB (+1) and any diagnostic heteroatom-specific neutral loss, HCN/CO/S/H2S
# (+2) - and drops the four weak/near-universal ones (KMD, M+1, M+2, aromatic
# DBE) entirely. Full per-criterion rationale is in R/30_wf3_hpac_filter.R.
wf3_scoring <- "minimal"     # "minimal", "recetox", or "original"

# Tier cut-offs are expressed as a FRACTION of the maximum achievable score so
# that all schemes use comparable tiers despite different maxima. Under the
# default "minimal" scheme (max 3): a diagnostic neutral loss alone (2/3 = 0.67)
# reaches High; a reference match alone (1/3 = 0.33) reaches Medium; both reach
# High - which matches the intended chemical priority (a heteroatom-specific
# neutral loss is stronger evidence than a bare fragment-mass match).
wf3_tier_high <- 0.55   # >= 55% of max  -> "High"
wf3_tier_med  <- 0.33   # >= 33% of max  -> "Medium"

# -----------------------------------------------------------------------------
# KMD model training source (R/50_kmd_confidence.R)
# -----------------------------------------------------------------------------
# "combined"        — fit the 12 KMD regression models from the empirical
#                      standards PLUS the full theoretical halo-PAC reference
#                      table (data_wf2/possible_halo_pacs2.csv). KM and KMD are
#                      deterministic functions of mass alone, so a theoretical
#                      entry contributes exactly the same information a
#                      physical standard of that structure would; see the full
#                      justification and the empirical checks that supported
#                      this choice in the header comment of
#                      load_theoretical_kmd_training() in that file.
# "standards_only"   — fit using only the physically-run standards, as before.
#                      Use this to reproduce earlier thesis/manuscript numbers
#                      exactly, or to compare against the combined result.
#
# Either setting only changes what trains the regression LINES. The points
# actually drawn on the figures (the grey training-standard triangles) are
# always the physical standards only, never the theoretical entries, so the
# plots are never cluttered with the ~100 additional theoretical points.
kmd_training_source <- "combined"   # "combined" or "standards_only"

# Exact monoisotopic neutral losses used as heteroatom diagnostics.
NEUTRAL_LOSSES <- c(
  HCN = 27.0109,   # hydrogen cyanide  - N-heterocycle diagnostic
  CO  = 27.9949,   # carbon monoxide   - O-heterocycle / carbonyl diagnostic
  S   = 31.9721,   # atomic sulfur     - thio-PAC diagnostic
  H2S = 33.9877,   # hydrogen sulfide  - thio-PAC diagnostic
  OH  = 17.0027    # hydroxyl radical  - phenolic / oxa-PAC diagnostic
)

# -----------------------------------------------------------------------------
# KMD series definitions (used by WF2 diagnostics and the final KMD stage)
# -----------------------------------------------------------------------------
# For a repeat unit with exact mass e and nominal mass n:
#   KM  = observed_mass * (n / e)
#   KMD = round(KM) - KM
# Members of the same homologous series share a KMD value.
#
# NOTE on naming: the "Cl" series uses the nominal/exact mass of a single
# chlorine atom (34 / 33.96103), which represents one H -> Cl substitution.
# It is labelled "H/Cl" rather than "Cl2" for that reason.
kmd_series <- list(
  CH2 = list(exact = 14.01565,  nominal = 14L,  label = "Kendrick Mass (CH\u2082 base)"),
  Cl  = list(exact = 33.96103,  nominal = 34L,  label = "Kendrick Mass (H/Cl base)"),
  Br2 = list(exact = 157.83604, nominal = 158L, label = "Kendrick Mass (Br\u2082 base)")
)

# Halogen classes that receive their own KMD regression model.
target_classes <- c("mono_Cl", "di_Cl", "mono_Br", "di_Br")

class_colours <- c(mono_Cl = "#1f77b4", di_Cl = "#5ba3d0",
                   mono_Br = "#d62728", di_Br = "#f4a582")

# -----------------------------------------------------------------------------
# Monoisotopic atomic masses (IUPAC / NIST), used to compute exact masses from
# molecular formulae during reference-library validation.
# -----------------------------------------------------------------------------
ATOMIC_MONO <- c(
  C = 12.0000000, H = 1.0078250, N = 14.0030740, O = 15.9949146,
  S = 31.9720707, Cl = 34.9688527, Br = 78.9183376, F = 18.9984032,
  I = 126.9044730, P = 30.9737615, Si = 27.9769265,
  D = 2.0141018   # deuterium, for labelled internal standards
)

message("Setup complete | sample = '", sample_name,
        "' | WF3 scoring = '", wf3_scoring, "'")
