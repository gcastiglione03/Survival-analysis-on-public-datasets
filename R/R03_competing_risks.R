# ==============================================================================
# 03_competing_risks.R
# Competing risks analysis: cumulative incidence function, Gray's test,
# Fine-Gray subdistribution hazard regression
# Author: Giuseppe Castiglione
# Dataset: MASS::Melanoma
# ==============================================================================

library(here)
library(dplyr)
library(ggplot2)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)

# DATA PREP --------------------------------------------------------------------

# status: 1=died from melanoma, 2=alive, 3=dead from other causes
# Recode to standard competing risks convention:
# 0=alive (censored), 1=melanoma death (event of interest), 2=other-cause death (competing event)
data(Melanoma, package = "MASS")

Melanoma <- Melanoma |>
  mutate(status = as.factor(case_when(
    status == 2 ~ 0,
    status == 1 ~ 1,
    status == 3 ~ 2
  )))

# WHY COMPETING RISKS ----------------------------------------------------------

# Applying 1-KM to melanoma-specific death treats other-cause deaths as censored,
# implicitly assuming they could still die from melanoma. This inflates the
# cumulative incidence. The CIF correctly accounts for the competing event.

# CUMULATIVE INCIDENCE FUNCTION ------------------------------------------------

# Both event types on one plot
cif_plot <- cuminc(Surv(time, status) ~ 1, data = Melanoma) |>
  ggcuminc(outcome = c("1", "2")) +
  ylim(c(0, 1)) +
  labs(
    x = "Days",
    title = "Cumulative incidence: melanoma death vs other-cause death"
  )

ggsave(here("outputs", "03_cif_overall.png"), cif_plot,
       width = 8, height = 6, dpi = 150)

# Output: the two curves sum to less than 1 at any timepoint, the remainder is
# the probability of still being alive. This is a key property of the CIF that
# 1-KM does not respect: applying 1-KM separately to each cause would produce
# curves that sum to more than 1, which is mathematically impossible and
# clinically misleading.


# CIF by ulceration status with Gray's test
# Gray's test: tests for difference in CIF between groups across full follow-up
cuminc(Surv(time, status) ~ ulcer, data = Melanoma) |>
  tbl_cuminc(
    times = 1826.25,
    label_header = "**{time/365.25}-year cuminc**"
  ) |>
  add_p()

cif_ulcer_plot <- cuminc(Surv(time, status) ~ ulcer, data = Melanoma) |>
  ggcuminc() +
  labs(
    x = "Days",
    title = "Cumulative incidence of melanoma death by ulceration status"
  ) +
  add_confidence_interval() +
  add_risktable()

ggsave(here("outputs", "03_cif_by_ulcer.png"), cif_ulcer_plot,
       width = 8, height = 6, dpi = 150)

# Output: Gray's test is the competing risks analogue of the log-rank test.
# It tests for differences in the CIF between groups across the full follow-up,
# accounting for the competing event. A significant p-value (p<0.001) here means
# the cumulative incidence of melanoma death differs by ulceration status, not
# just the hazard rate.


# FINE-GRAY REGRESSION ---------------------------------------------------------

# crr() models the subdistribution hazard for the event of interest (among all patients)
# Subjects who experience the competing event remain in the risk set (with downweighting)
crr(Surv(time, status) ~ sex + age, data = Melanoma) |>
  tbl_regression(exp = TRUE)

# Alternative: cause-specific Cox model (censors competing events, melanoma death only)
# Use when the clinical question is about the hazard mechanism, not absolute risk
coxph(
  Surv(time, ifelse(status == 1, 1, 0)) ~ sex + age,
  data = Melanoma
) |>
  tbl_regression(exp = TRUE)