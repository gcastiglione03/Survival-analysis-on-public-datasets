# ==============================================================================
# 01_kaplan_meier.R
# Kaplan-Meier estimation, log-rank test, Cox proportional hazards regression
# Author: Giuseppe Castiglione
# Dataset: survival::lung
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(gtsummary)


# DATA PREP --------------------------------------------------------------------

# lung dataset: time (days), status (1=censored, 2=dead), sex (1=Male, 2=Female)
# Recode to standard convention: 0=censored, 1=dead
lung <- lung |>
  mutate(status = case_when(
    status == 1 ~ 0,
    status == 2 ~ 1
  ))

# Format dates for R
date_ex <-
  tibble(
    sx_date = c("2007-06-22", "2004-02-13", "2010-10-27"),
    last_fup_date = c("2017-04-15", "2018-07-04", "2016-10-31")
  )
date_ex


# KAPLAN-MEIER CURVE -----------------------------------------------------------

# Overall survival curve (no grouping)
km_fit <- survfit2(Surv(time, status) ~ 1, data = lung)

km_plot <- km_fit |>
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability",
    title = "Overall survival — NCCTG Lung Cancer Trial"
  ) +
  add_confidence_interval() +
  add_risktable()

ggsave(here("outputs", "01_km_overall.png"), km_plot,
       width = 8, height = 6, dpi = 150)

# Output: the curve drops steadily with no plateau (especially in the first 200 days), suggesting there is no
# long-term survivor subpopulation in this cohort.


# 1-YEAR SURVIVAL --------------------------------------------------------------

# Correct estimate: 40.9% (naive estimate ignoring censoring before 1yr would give ~47%)
survfit(Surv(time, status) ~ 1, data = lung) |>
  tbl_survfit(
    times = 365.25,
    label_header = "**1-year survival (95% CI)**"
  )

# Output: fewer than half of patients survive to one year.
# A naive estimate (deaths / total) would yield ~47% because it ignores
# patients censored before 365 days. KM corrects for this by treating
# censoring as informative about time, not about outcome.


# MEDIAN SURVIVAL --------------------------------------------------------------

# Median survival: 310 days (50% of patients have already died)
survfit(Surv(time, status) ~ 1, data = lung) |>
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )

# Median is the preferred summary statistic in oncology over the mean, because
# survival time distributions are typically right-skewed (a few long survivors
# would inflate the mean artificially).


# SURVIVAL BY SEX --------------------------------------------------------------

km_sex_plot <- survfit2(Surv(time, status) ~ sex, data = lung) |>
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability",
    title = "Survival by sex — NCCTG Lung Cancer Trial"
  ) +
  add_confidence_interval() +
  add_risktable()

ggsave(here("outputs", "01_km_by_sex.png"), km_sex_plot,
       width = 8, height = 6, dpi = 150)

# Log-rank test: female patients show significantly longer survival (p < 0.001)
survdiff(Surv(time, status) ~ sex, data = lung)

# Output: statistically significant (p < 0.001) difference between the two groups.
# The log-rank test compares curves across the entire follow-up period,
# not at a single timepoint (better than a 1-year comparison because it uses all
# available time-to-event information.


# COX REGRESSION ---------------------------------------------------------------

# HR < 1 for sex confirms survival advantage for women (sex=2)
coxph(Surv(time, status) ~ sex, data = lung) |>
  tbl_regression(exp = TRUE)

# Output: HR = 0.59 for sex. Women (sex=2) have a 41% lower risk of death
# compared to men, at any given point in follow-up.
# The Cox model quantifies this effect as a hazard ratio assumed constant
# over time.
