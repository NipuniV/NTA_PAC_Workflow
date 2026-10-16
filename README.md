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
