# Survival Analysis in R
**Time-to-event analysis from first principles to competing risks** — Kaplan-Meier estimation, landmark analysis, time-dependent covariates, Fine-Gray subdistribution hazards, and conditional survival, implemented across four progressively complex clinical datasets.

---

## Overview

This project implements a complete survival analysis workflow in R, starting from basic Kaplan-Meier curves and building toward competing risks regression and conditional survival estimation. Each section addresses a specific methodological challenge that arises in real clinical data — not just how to run the code, but when and why each approach is appropriate.

The datasets used are publicly available: `lung` (NCCTG Lung Cancer Trial), `BMT` (bone marrow transplant, 137 patients), and `Melanoma` (cutaneous melanoma, {MASS} package).

---

## Methods

| Section | Dataset | Method | Key question addressed |
|---------|---------|--------|------------------------|
| Part 1 | `lung` | Kaplan-Meier, log-rank test, Cox PH | Estimating and comparing survival curves |
| Part 2 | `BMT` | Landmark analysis, time-dependent covariates | Handling post-baseline events (acute GVHD) |
| Part 3 | `Melanoma` | Competing risks (Fine-Gray), Gray's test | When 1 − KM overestimates the cumulative incidence |
| Part 4 | `lung` | PH assumption testing, conditional survival | Model diagnostics and survivorship bias |

---

## Why competing risks matter

The standard Kaplan-Meier estimator treats all non-events as censored — including deaths from competing causes. In the Melanoma dataset, patients can die from melanoma (event of interest) or from other causes (competing event). Applying 1 − KM to the melanoma-specific endpoint inflates the cumulative incidence because it implicitly assumes the competing risk never occurs.

The correct approach — implemented here via `tidycmprsk::cuminc()` and `crr()` — uses the **cumulative incidence function (CIF)**, which properly accounts for the competing cause. The difference is not academic: in the 5-year melanoma data, the two methods produce estimates that diverge meaningfully enough to affect clinical interpretation.

---

## Repository structure

```
survival-analysis-r/
├── R/
│   ├── 01_kaplan_meier.R          # KM curves, log-rank, Cox regression
│   ├── 02_landmark_time_dependent.R   # Landmark analysis, tmerge, TD covariates
│   ├── 03_competing_risks.R       # CIF, Gray's test, Fine-Gray regression
│   └── 04_advanced_topics.R       # PH diagnostics, smooth survival, conditional survival
├── outputs/                       # Saved plots (.png)
├── .gitignore
└── README.md
```

---

## Setup

All packages are installed from CRAN except `condsurv`, which is installed from GitHub.

```r
install.packages(c(
  "dplyr", "ggplot2", "lubridate",
  "survival", "ggsurvfit", "gtsummary",
  "tidycmprsk", "SemiCompRisks", "MASS",
  "sm", "tibble"
))

pak::pak("zabore/condsurv")
```

R version used: **4.3.x**. Package versions are tracked in `renv.lock` (see below).

---

## Key outputs

### Kaplan-Meier with risk table — `lung` dataset
- 1-year overall survival: **40.9%** (95% CI: 34–48%)
- Median survival: **310 days**
- Female patients showed significantly longer survival (log-rank p < 0.001, HR = 0.59)

### Competing risks — `Melanoma` dataset
- 5-year cumulative incidence of melanoma-specific death: differentiated by ulceration status (Gray's test p < 0.05)
- Fine-Gray model: sex and age as predictors of the subdistribution hazard

### Conditional survival — `lung` dataset
- Among patients alive at 6 months, estimated probability of surviving to 12 months
- Implemented via `condsurv::conditional_surv_est()` across multiple conditional time points

---

## Methodological notes

**Landmark analysis vs time-dependent covariates.** Landmark analysis is simpler but requires a pre-specified cutoff and excludes patients who died before it (15 patients excluded at 90 days in the BMT data). Time-dependent covariates via `tmerge()` avoid this selection bias and are preferred when the landmark is arbitrary or exclusions are substantial.

**Proportional hazards assumption.** Tested with `cox.zph()` (Schoenfeld residuals). Both `sex` and `age` satisfied the PH assumption in the lung dataset (p > 0.05), validating the use of a standard Cox model.

**Censoring.** Status variables were recoded to the standard convention (0 = censored, 1 = event) at the start of each analysis. In Part 3, status = 0 (alive), 1 (melanoma death), 2 (other-cause death) — the competing event.

---

## Skills demonstrated

`R` · `survival` · `ggsurvfit` · `tidycmprsk` · `gtsummary` · Kaplan-Meier estimation · Cox proportional hazards regression · Fine-Gray subdistribution hazard model · competing risks analysis · cumulative incidence function (CIF) · landmark analysis · time-dependent covariates · `tmerge` · conditional survival · Gray's test · Schoenfeld residuals · time-to-event analysis · clinical biostatistics

---

## References

Datasets: `survival::lung`, `SemiCompRisks::BMT`, `MASS::Melanoma`, `condsurv` package (Zabore, GitHub).

---

## Author

**Giuseppe Castiglione** — Bioinformatician & Data Scientist, Brussels  
[LinkedIn](https://linkedin.com/in/giuseppecastiglione03) · [GitHub](https://github.com/gcastiglione03)
