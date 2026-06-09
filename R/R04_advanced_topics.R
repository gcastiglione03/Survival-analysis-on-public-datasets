# ==============================================================================
# 04_advanced_topics.R
# Advanced topics: proportional hazards assumption, smooth survival curves,
# conditional survival estimation
# Author: Giuseppe Castiglione
# Dataset: survival::lung
# ==============================================================================

library(here)
library(dplyr)
library(survival)
library(ggsurvfit)
library(sm)

pak::pak("zabore/condsurv")
library(condsurv)

# DATA PREP --------------------------------------------------------------------

lung <- lung |>
  mutate(status = case_when(
    status == 1 ~ 0,
    status == 2 ~ 1
  ))

# PROPORTIONAL HAZARDS ASSUMPTION ----------------------------------------------

# The Cox model assumes hazards are proportional at every timepoint.
# cox.zph() tests this via Schoenfeld residuals.
# p > 0.05 = do not reject PH assumption for that covariate.
mv_fit <- coxph(Surv(time, status) ~ sex + age, data = lung)
cz <- cox.zph(mv_fit)
print(cz)
# Output: this validates the Cox models fitted in 01_kaplan_meier.R for both sex 
# and age in this dataset. The HR of 0.59 for sex can be interpreted as constant
# across the entire follow-up period, not just at a single timepoint.

ph_plot_path <- here("outputs", "04_ph_assumption.png")
png(ph_plot_path, width = 800, height = 500, res = 120)
plot(cz)
dev.off()


# SMOOTH SURVIVAL (CONTINUOUS COVARIATE) ---------------------------------------

# sm.survival() estimates median survival time as a smooth function of age.
# Bandwidth h chosen by Silverman's rule-of-thumb.
smooth_plot_path <- here("outputs", "04_smooth_survival_age.png")
png(smooth_plot_path, width = 800, height = 500, res = 120)

sm.options(list(
  xlab = "Age (years)",
  ylab = "Median time to death (years)"
))

sm.survival(
  x      = lung$age,
  y      = lung$time,
  status = lung$status,
  h      = (1/6) * sd(lung$age) / nrow(lung)^(-1/4)
)

dev.off()

# Output: # The smooth curve reveals that the relationship between age and
# median survival time is linear. Median survival time decreases when age increases.
#Age only is a weak predictor of survival in this dataset.


# CONDITIONAL SURVIVAL ---------------------------------------------------------

# Conditional survival: P(survive t2 | already survived t1)
# Addresses survivorship bias, patients who have already survived some time
# face a different risk profile than the original cohort.

km_fit <- survfit(Surv(time, status) ~ 1, data = lung)

# Estimate conditional survival at 6-month intervals up to 2 years,
# conditional on having survived the first 6 months (t1 = 182.625 days)
prob_times <- seq(365.25, 182.625 * 4, 182.625)

conditional_estimates <- purrr::map_df(
  prob_times,
  ~conditional_surv_est(
    basekm = km_fit,
    t1     = 182.625,
    t2     = .x
  )
) |>
  mutate(months = round(prob_times / 30.4)) |>
  select(months, everything())

print(conditional_estimates)

# Visualization across multiple conditional time points
cond_plot <- gg_conditional_surv(
  basekm = km_fit,
  at     = prob_times,
  main   = "Conditional survival — NCCTG Lung Cancer Trial",
  xlab   = "Days"
) +
  labs(color = "Conditional time (days)")

ggsave(here("outputs", "04_conditional_survival.png"), cond_plot,
       width = 8, height = 6, dpi = 150)

# Output: the plot shows three conditioning times (1yr, 18mo, 2yr).
# Curves shift upward with increasing conditioning time: patients who have
# already survived longer show progressively better subsequent survival,
# the quantitative signature of survivorship bias.
# The blue curve (t1=730 days) stabilizes around 45% but is based on few
# patients still at risk beyond day 1000, so estimates carry high uncertainty.
# This pattern is clinically relevant for patient counseling at follow-up:
# the unconditional KM curve from 01_kaplan_meier.R (median = 310 days)
# substantially underestimates survival prospects for patients already
# alive at 1 or 2 years post-diagnosis.