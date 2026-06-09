# ==============================================================================
# 02_landmark_time_dependent.R
# Landmark analysis and time-dependent covariates
# Author: Giuseppe Castiglione
# Dataset: SemiCompRisks::BMT (137 bone marrow transplant patients)
# ==============================================================================

library(here)
library(dplyr)
library(tibble)
library(survival)
library(ggsurvfit)
library(gtsummary)

# DATA PREP --------------------------------------------------------------------

# T1: time to death or last follow-up (days)
# delta1: death indicator (1=Dead, 0=Alive)
# TA: time to acute graft-versus-host disease (aGVHD)
# deltaA: aGVHD indicator (1=Developed, 0=Never developed)
data(BMT, package = "SemiCompRisks")
BMT <- rowid_to_column(BMT, "my_id")

# LANDMARK ANALYSIS ------------------------------------------------------------

# Landmark time: 90 days (chosen on clinical grounds, prior to data inspection)
# 15 patients excluded: died before reaching the 90-day landmark
lm_dat <- BMT |>
  filter(T1 >= 90) |>
  mutate(lm_T1 = T1 - 90)   # follow-up time measured from landmark

# NOTE: this is the core limitation of landmark analysis: it conditions the
# entire analysis on surviving to the landmark, introducing selection bias if
# the excluded patients differ systematically from those who remain.# 

lm_plot <- survfit2(Surv(lm_T1, delta1) ~ deltaA, data = lm_dat) |>
  ggsurvfit() +
  labs(
    x = "Days from 90-day landmark",
    y = "Overall survival probability",
    title = "Survival from 90-day landmark by aGVHD status"
  ) +
  add_risktable()

ggsave(here("outputs", "02_landmark_plot.png"), lm_plot,
       width = 8, height = 6, dpi = 150)

# Output: patients who developed aGVHD before day 90 show a different survival
# trajectory from those who did not. However, this comparison is already
# conditioned on surviving to the landmark. The effect size may be
# underestimated if early aGVHD deaths were excluded disproportionately.


# Cox model restricted to patients surviving past landmark
coxph(Surv(T1, delta1) ~ deltaA, subset = T1 >= 90, data = BMT) |>
  tbl_regression(exp = TRUE)

# Output: the subset argument replicates the landmark restriction inside coxph()
# without creating a new dataset. The HR (1.08) quantifies the association between
# aGVHD and mortality among patients who survived at least 90 days, not in the 
# full original cohort.


# TIME-DEPENDENT COVARIATE APPROACH --------------------------------------------

# Preferred over landmark when:
# - no obvious landmark time exists
# - landmark exclusions are substantial (here: 15/137 patients)
# - the covariate value changes over time

# tmerge() creates a long-format dataset with one row per time interval per patient
# event() defines the outcome within each interval
# tdc() defines the time-dependent covariate (aGVHD onset)
td_dat <- tmerge(
  data1 = BMT |> select(my_id, T1, delta1),
  data2 = BMT |> select(my_id, T1, delta1, TA, deltaA),
  id = my_id,
  death = event(T1, delta1),
  agvhd = tdc(TA)
)

# Surv() with time + time2 arguments handles the interval-censored structure
coxph(
  Surv(time = tstart, time2 = tstop, event = death) ~ agvhd,
  data = td_dat
) |>
  tbl_regression(exp = TRUE)

# NOTE: the two-argument Surv(tstart, tstop, event) form tells Cox to evaluate
# the hazard only within each interval, where agvhd has a fixed value.
# Compared to the landmark model, this uses all 137 patients and avoids
# the 90-day selection cutoff. The HR estimate is therefore less biased
# and based on more complete information.