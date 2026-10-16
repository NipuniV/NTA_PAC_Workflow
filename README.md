# Integrated Non-Targeted Analysis Workflow for Polycyclic Aromatic Compounds

A single, automated pipeline for identifying homocyclic, heterocyclic and
halogenated PACs in GCxGC-HRTOF-MS data. It consolidates five previously
separate R projects into one reproducible workflow that runs end to end from a
single command, while keeping Workflows 1-4 clearly and separately defined.

---

## Quick start

1. Place your ChromaTOF export (e.g. `YP1.csv`) into **`data_raw/`**.
   It goes there **once** — Workflows 1, 2 and 3 all read the same file.
2. Open `NTA_PAC_Workflow.Rproj` in RStudio.
3. Set the sample name near the top of `Rmarkdown/00_setup.Rmd`:
   ```r
   sample_name <- "YP1"
   ```
4. Knit `run_all.Rmd` (the **Knit** button, or
   `rmarkdown::render("run_all.Rmd")`).

Knitting produces `run_all.html`: a single report showing what every stage did,
how many peaks survived each filter, and where each file was written.

Any required package that is not already installed is installed automatically by
the setup stage, so the project runs on a clean machine with no preparation.

Primary outputs:

| File | Contents |
|---|---|
| `results_final/<sample>_final_PAC_report.csv` | Integrated WF4 report |
| `results_final/<sample>_high_confidence_halo_PACs.csv` | Tiered WF2 + KMD halo-PAC list |
| `kmd_final/standards_kmd_HCl.png` | H/Cl models, training data only |
| `kmd_final/standards_kmd_Br2.png` | Br2 models, training data only |
| `kmd_final/standards_kmd_CH2.png` | CH2 models, training data only |
| `kmd_final/<sample>_kmd_HCl.png` | H/Cl series KMD plot |
| `kmd_final/<sample>_kmd_Br2.png` | Br2 series KMD plot |
| `kmd_final/<sample>_kmd_CH2.png` | CH2 series KMD plot |

---

## Folder layout

```
NTA_PAC_Workflow/
├── run_all.Rmd                     knit this — runs every stage in order
├── Rmarkdown/
│   ├── 00_setup.Rmd                ALL tolerances and thresholds live here
│   ├── 01_helpers.Rmd              functions shared between workflows
│   ├── 02_validate_libraries.Rmd   automated reference-library integrity check
│   ├── 10_wf1_library_match.Rmd    WORKFLOW 1
│   ├── 20_wf2_halogen_filter.Rmd   WORKFLOW 2
│   ├── 30_wf3_hpac_filter.Rmd      WORKFLOW 3
│   ├── 40_wf4_integration.Rmd      WORKFLOW 4
│   └── 50_kmd_confidence.Rmd       KMD models + final plots
├── R/                              plain-script copies of the same code
├── data_raw/   the ChromaTOF export — stored ONCE, read by WF1, WF2 and WF3
├── data_wf1/   PAC_mass_spectral_library.csv
├── data_wf2/   possible_halo_pacs2.csv + std_check_all_halo_pacs.csv
├── data_wf3/   NSO_PAC_exact_mass_reference.csv + HPAC_mass_spectral_library.csv
├── data_wf4/   outputs of WF1, WF2 and WF3 (staged automatically)
├── results_wf1..4/   every intermediate file each workflow writes
├── results_final/    the integrated WF4 report
├── kmd_final/        the three combined KMD plots
└── logs/             library validation reports
```

**Sample data lives in one place.** `data_raw/` holds the ChromaTOF export, and
Workflows 1, 2 and 3 all read it from there — there is no duplicate copy to keep
in sync. The `data_wfN/` folders hold only the reference tables that workflow
needs.

Each workflow writes only to its own `results_wfN/` folder. Files needed by a
later stage are copied into `data_wf4/` automatically.

**Every workflow is a standalone document.** Knitting
`Rmarkdown/20_wf2_halogen_filter.Rmd` on its own produces a self-contained
report for Workflow 2: it loads the configuration and helpers it needs, runs,
and reports. When knitted as part of `run_all.Rmd` it skips its own run chunk so
nothing executes twice. The `R/` folder keeps plain-script copies of the same
code for anyone who prefers `source()` over knitting.

---

## The four workflows

WF1, WF2 and WF3 are **independent** and interrogate the same peak list from
three different angles. Only WF4 depends on the others.

| | Question asked | Nature of the evidence |
|---|---|---|
| **WF1** | Does this peak match an authentic standard? | Spectral + retention match |
| **WF2** | Does this peak contain Cl or Br? | Isotope ratios — a physical constant |
| **WF3** | Does this peak contain a ring N, S or O? | Diagnostic neutral losses, mass defect |
| **WF4** | What do all three together support? | Convergence of independent evidence |

### Why WF2 is strong and WF3 is deliberately not

The isotope abundances of `37Cl` and `81Br` are fixed properties of the nuclide,
so a well-formed isotope cluster is close to decisive evidence of halogen
content. Nitrogen, sulfur and oxygen offer no equivalent signature under EI:
their mass-defect contributions are small, and the neutral losses they produce
(HCN, CO, H2S) are shared with many ordinary organics.

WF3 is therefore built for **high recall, not high precision**. An HPAC dropped
at the screening stage can never be recovered; a false candidate can be removed
later by WF4. Its output is a candidate list, never a set of identifications.

---

## Two things worth knowing before you run it

### 1. Reference libraries are validated automatically, before any data is processed

Every reference table is checked by recomputing each compound's monoisotopic
mass from its own molecular formula and comparing it to the tabulated value.

This exists because of a real failure. The reference entry for
1,6-dibromopyrene carried **360.0419 Da** — the *average molecular weight* of
C16H8Br2, not its monoisotopic mass (**357.8993 Da**). WF2 correctly recognised
the compound as di-brominated from its isotope pattern, but the exact-mass
matching step could never assign it, because the reference value sat 2.14 Da
away from anything the instrument could measure. There was no error message:
the compound simply vanished from the results and had to be confirmed by hand.

Running the same check across all four reference tables found **31 further
entries** with the same class of problem. Corrected libraries are included, and
every change is logged in `logs/reference_library_corrections.csv` with the
original value, the corrected value and the reason.

Two entries needed judgement rather than arithmetic and are flagged in that log:

- `benzo[b]fluoranthene` carried the formula `C20D12` (the deuterated
  standard's formula, copy-pasted). Corrected to `C20H12`.
- `benzo[k]fluoranthene-d13` is named d13 but carries `C20D12`. Benzo[k]-
  fluoranthene is C20H12, so d12 is the maximum labelling and the name is most
  likely a typo. **The formula was left unchanged — please confirm.**

`possible_halo_pacs2.csv` (300 entries) and `NSO_PAC_exact_mass_reference.csv`
(154 entries) passed with no errors and are unchanged.

### 2. KMD training standards use verified ground truth, not WF2's own guess

The KMD models (Kendrick mass defect regression lines used to score unknown
compounds) are trained on the 28 authentic standards in the mixture. Their
identity is already known with certainty — they are not unknowns — so the
KMD stage does **not** use WF2's own isotope-pattern classification to decide
which training group a standard belongs to. It instead reads
`data_wf2/std_check_ground_truth_halo_class.csv`, a copy of WF2's output in
which `halogen_class` has been manually verified against each compound's
known molecular formula.

This exists because WF2's own classification is not always right, and using
a wrong guess to label a *training* example is worse than using no label at
all — it teaches the model the wrong pattern. In this project's own
standards, WF2 placed 1,6-Br2-Pyr in the mono_Cl group rather than di_Br,
purely because of which isotope peaks it happened to pair together, despite
the compound's identity never being in doubt. Several other dibrominated
standards were affected the same way.

Real, unknown sample compounds are unaffected by this: for those, WF2's
isotope-pattern evidence remains the only tool available and is used exactly
as before.

**This ground-truth file is itself checked automatically**, before every use,
by `validate_standards_ground_truth()` in `R/02_validate_libraries.R`. It
recomputes each standard's expected class directly from its known formula
(independent of WF2 entirely) and flags any manually-entered class that
disagrees. This check is what caught `5,6-Br2-1,2-H2-Acy` entered as
`poly_or_mixed` when its formula (C12H8Br2) is unambiguously `di_Br` during
preparation of this file — fix any such entries the check reports before
trusting the KMD models built from it.

### 3. KMD regression models train on more than just the physical standards

The 12 KMD regression models (one per series × halogen class) can be fitted
from the physical standards alone, or from those standards **combined with**
the full theoretical halo-PAC reference table
(`data_wf2/possible_halo_pacs2.csv`). Controlled by `kmd_training_source` in
`R/00_setup.R`:

```r
kmd_training_source <- "combined"   # or "standards_only"
```

**Why theoretical entries are valid training data at all.** Kendrick mass and
Kendrick mass defect are both deterministic functions of mass alone
(`KM = mass * nominal/exact`, `KMD = round(KM) - KM`). Neither depends on how
that mass was obtained. A regression of "where a halogen class sits in KMD
space" is a model of chemistry — how the H:C ratio trends with molecular size
within a class — not of measurement, so there is no principled reason to
restrict training data to the handful of physical standards per class when
up to 25 theoretical examples per class exist in the same reference table WF2
already uses for exact-mass matching.

**This was checked empirically, not assumed.** Fitting each class from
theory-only and standards-only data separately gave closely matching slopes
and comparable residual scatter for three of the four classes (mono_Cl,
di_Cl, mono_Br) — real measurement noise, which theoretical masses cannot
capture, was checked for directly and found not to be the dominant source of
scatter in any of them. The fourth class (di_Br) showed an apparent slope
disagreement that turned out to be small-sample instability rather than a
real conflict: the theoretical slope fell well inside the standards-only
slope's own 95% confidence interval, and combining the two reduced the
standard error of the fitted slope by 85%.

| Class | n, standards only | n, combined |
|---|---|---|
| mono_Cl | 5 | 30 |
| di_Cl | 4 | 29 |
| mono_Br | 10 | 35 |
| di_Br | 5 | 30 |

**Classification is by formula, never by name.** An earlier attempt
classified theoretical entries by checking whether the compound name
contained the text "Mono-Cl" or "Di-Br" etc. This silently mis-classified
mixed-halogen entries such as `Indole-Mono-Cl-Di-Br` as `mono_Cl`, because
that substring happens to appear in the name despite the compound containing
both chlorine and bromine. `load_theoretical_kmd_training()` in
`R/50_kmd_confidence.R` instead counts Cl and Br atoms directly from the
molecular formula, which cannot make this mistake.

**The figures are unaffected by this choice.** Only the fitted regression
lines and their prediction-interval bands change; the points actually drawn
on each plot (the grey training-standard triangles) are always the physical
standards only. The ~100 additional theoretical points inform the model
without appearing individually in the figures.

Use `"standards_only"` to reproduce numbers from before this feature existed.

**Model fitting is now a separate step from testing a sample**, mirroring the
original standalone KMD project this pipeline consolidates. `build_kmd_models()`
fits the 12 regression models and immediately saves a standards-only
diagnostic figure per series — `kmd_final/standards_kmd_HCl.png`,
`..._Br2.png`, `..._CH2.png` — showing the fitted lines and the physical
training standards with **no sample data at all**. This runs once, before any
sample is processed, so the models can be inspected on their own. `run_kmd(sample_name)`
then tests that specific sample's WF2 output against the already-built models
and saves `kmd_final/<sample>_kmd_*.png` as before.

`run_all.R` and `run_all.Rmd` call both steps explicitly and in order. Calling
`run_kmd(sample_name)` on its own (e.g. when knitting `50_kmd_confidence.Rmd`
standalone) still works exactly as before — it builds the models itself first
if none are supplied.

### 4. ChromaTOF exports with unlabelled columns are handled

Some exports carry trailing columns that have **no header label but do contain
ion data** — YP1 has 72 such columns holding data for 752 spectra, and 9031_B1
has 122 holding data for 381. They are the 19th ion onwards for spectra with
more ions than the header had labels for.

These files also frequently contain Excel error strings (`#VALUE!`, `#SPILL!`)
left in the ion columns.

WF2 now reads the export entirely as text and coerces every column explicitly,
rather than letting `readr` infer types from the first 1000 rows. Ion columns
are paired by number where they are labelled, and the unlabelled tail is paired
positionally so none of that data is lost. Verified to reproduce the original
per-sample results exactly:

| Sample | Unlabelled columns | Peaks passing the halogen filter |
|---|---|---|
| `std_check` | 0 | 28 |
| `std1_LE` | 0 | 92 |
| `YP1` | 72 | 527 |
| `9031_B1` | 122 | 70 |

### 5. WF3 has three scoring schemes — pick deliberately

Set in `R/00_setup.R`:

```r
wf3_scoring <- "minimal"   # "minimal" (default), "recetox", or "original"
```

**`"minimal"` (default, manuscript stage)** keeps only the two criteria with
genuine, defensible chemical discriminatory power and drops the four weak ones
entirely rather than down-weighting them. Max 3 points:

| Criterion | Points | Why it is kept |
|---|---|---|
| Fragment-ion match to PAC/HPAC reference DB (±5 ppm) | +1 | Under EI, PACs fragment before the molecular ion is recorded, so matching *all* detected fragment masses (not just the base peak) to known PAC/HPAC references still catches diagnostic substructure fragments. Anchors a candidate to PAC-like fragmentation. Moderate discrimination. |
| Any diagnostic neutral loss (HCN / CO / S / H₂S) | +2 | Heteroatom rings fragment along heteroatom-specific pathways — N→HCN, O→CO, S→S/H₂S — so the lost neutral maps directly onto which heteroatom is present. The only criterion shown to meaningfully separate HPACs from background. Strong discrimination. |

Dropped: `kmd_flag` (fired on 93.9% of all features — flags "alkylated" more
than "heteroatom-containing"), `has_M1` (present in essentially every
carbon-containing compound), `has_M2` (genuine S pointer but shared with any
S-compound; slightly negative discrimination), `dbe_aromatic` (fired on 99.2%
of features — true of homocyclic PAHs too). Tier mapping under this scheme: a
diagnostic neutral loss alone (2/3) reaches **High**; a reference match alone
(1/3) reaches **Medium**; both reach **High** — matching the intended priority.

**`"recetox"`** applies a re-weighting derived from benchmarking WF3 against the
RECETOX Exposome library (386 spectra, 23 true HPACs). Each criterion was
measured for *discrimination* — how much more often it fires for a true HPAC
than for background:

| Criterion | Fires overall | Discrimination |
|---|---|---|
| `kmd_flag` (2 points) | 93.9% | +0.018 |
| `dbe_aromatic` | 99.2% | +0.009 |
| `has_M1` | 93.4% | +0.071 |
| `has_M2` | 16.0% | **−0.031** |
| `sulfur_iso_flag` | 1.9% | **−0.021** |
| any neutral loss | 79.9% | +0.122 |
| `frag_N_loss` (HCN) | 35.0% | +0.137 |

It makes `dbe_aromatic` a mandatory gate, drops `has_M2`, and moves weight onto
the element-specific neutral losses. This was the intermediate step that
motivated the simpler `"minimal"` scheme now used by default.

**`"original"`** reproduces the thesis and submitted manuscript exactly (max 8
points). Use this whenever you need to regenerate published numbers.

**Neither `"minimal"` nor `"recetox"` reproduces published numbers.**

---

## The high-confidence halo-PAC list (`results_final/<sample>_high_confidence_halo_PACs.csv`)

A final integration step (`R/60_high_confidence_halo_pacs.R`, run after the KMD
stage; WF4 is untouched) combines the two halogen evidence streams into one
tiered list:

| Tier | What was independently established |
|---|---|
| **High** | WF2 isotope-positive **and** exact-mass reference match **and** KMD in-band for that class |
| **Medium (a)** | isotope + reference match, but KMD out-of-band → *class anomaly, review* |
| **Medium (b)** | isotope + KMD in-band, but no reference match → *possible novel halo-PAC* |
| **Low / review** | isotope pattern only |

**Why KMD is used for corroboration, not as independent proof.** The KMD models
are trained on the theoretical halo-PAC library, and Kendrick mass defect is a
function of mass alone with no ring-structure information. "Falls in the di-Br
band" therefore means "has a mass defect consistent with some di-brominated
composition" — which a di-brominated non-PAC could also satisfy. Treating
KMD-in-band as independent confirmation of *halo-PAC identity* would be partly
circular (rewarding a compound for resembling the set used to define the
classes) and overlaps the "is this halogenated?" signal WF2 already provides. So
this step uses KMD only for what it can establish without circularity: **class
corroboration** (does the KMD-supported class agree with WF2's?) and **anomaly
flagging** (a halogenated peak outside all bands is flagged for review, not
rejected). A `review_flag` column additionally marks High-tier peaks where KMD
supports a *different* class than WF2, and any KMD extrapolation warnings.

If the KMD stage produced no output, the step degrades gracefully: no peak can
reach High (which requires the KMD line), reference-matched peaks land in
Medium, and isotope-only peaks in Low.

## What changed from the original five projects

| Change | Reason |
|---|---|
| One command runs all stages | Previously each workflow was run by hand and its output copied into the next project |
| All thresholds in `00_setup.R` | The same tolerance previously appeared in several files and could drift out of sync |
| Shared functions in `01_helpers.R` | Each project carried its own copy of `add_kmd()`, the reshape logic and the AI calculation |
| Automated library validation | Catches the silent-failure class of error described above |
| Sample name set once | Was hardcoded separately in every script, with a comment to change it each run |
| Raw file stored once, in `data_raw/` | The same export previously had to be copied into three folders and kept in sync |
| R Markdown throughout | Running the pipeline now produces a readable HTML report, not just a console log |
| Packages installed automatically | The project runs on a clean machine without the user preparing anything first |
| Per-workflow data/results folders | Makes it obvious what each stage consumed and produced |
| WF3 stops after CSV export | The KMD and contour plots that previously followed are not produced; consolidated KMD analysis now lives in stage 5 |
| 15 KMD plots reduced to 3 | Only the combined plot per series is produced; per-class plots are no longer generated |

Every file the original scripts wrote is still written, so nothing that was
previously available for inspection has been lost.

---

## Running a single workflow

Either knit the workflow's own document —

```r
rmarkdown::render("Rmarkdown/20_wf2_halogen_filter.Rmd")
```

— or use the plain scripts:

```r
source("R/00_setup.R")
source("R/01_helpers.R")
source("R/20_wf2_halogen_filter.R")
run_wf2(sample_name)
```

Each workflow is a self-contained function (`run_wf1()` … `run_kmd()`) and can
be run, described and cited independently.

---

## Requirements

R (>= 4.1). Required packages are installed automatically on first run:
`tidyverse`, `here`, `janitor`, `ggrepel`, `knitr`, `rmarkdown`.

To install them yourself in advance:

```r
install.packages(c("tidyverse", "here", "janitor", "ggrepel",
                   "knitr", "rmarkdown"))
```
