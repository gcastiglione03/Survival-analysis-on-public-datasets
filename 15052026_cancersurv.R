################################################################################
#                                                                              #
# Project name: Survival analysis in R                                         #
# User: Giuseppe Castiglione                                                   #
# Creation date: 09/04/2026                                                    #
# Last edit date: 09/06/2026                                                   #
# Dataset: zabore/condsurv, package condsurv                                   #
#                                                                              #
################################################################################

# ==================
# Start
# ==================

# PACKAGES #--------------------------------------------------------------------

install.packages("here")
library(here)

library(dplyr)
library(ggplot2)
library(lubridate)
library(survival)
library(ggsurvfit)
library(gtsummary)
library(tidycmprsk)
library(tibble)

# install.packages("sm")
library(sm)

pak::pak("zabore/condsurv")
library(condsurv)


# DIRECTORY---------------------------------------------------------------------

PATH <- "C:/Users/..."
#NAME_INPUT <- "condsurv"
PATH_OUT <- paste0(PATH,"Survival/")
PATH_SCRIPT <- "C:/Users/..."

# PART 1: INTRODUCTION TO SURVIVAL ANALYSIS #-----------------------------------

# THE LUNG DATASET: format and date #-------------------------------------------
# time: Observed survival time in days
# status: censoring status 1=censored, 2=dead
# sex: 1=Male, 2=Female

# Recode the status in a standard way
lung <-
  lung |>
  mutate(
    status = case_when(status, `1` = 0, `2` = 1)
  )
# new status: censoring status 0=censored, 1=dead

head(lung[, c("time", "status", "sex")])

## Calculating survival times ##
# Format dates for R
date_ex <-
  tibble(
    sx_date = c("2007-06-22", "2004-02-13", "2010-10-27"),
    last_fup_date = c("2017-04-15", "2018-07-04", "2016-10-31")
    )
date_ex


# SURVIVAL OBJECTS AND CURVES #-------------------------------------------------

# Creating survival objects via Surv()
Surv(lung$time, lung$status)[1:10]
# Output: survival time, followed by a "+" if the subject was censored

# survfit(): creates survival curves using the Kaplan-Meier method
# No categories to compare so in the function of nothing, so "~ 1"
s1 <- survfit(Surv(time, status) ~ 1, data = lung)
str(s1)
# Output
# time: the timepoints at which the curve has a step, i.e. at least one event occurred
# surv: the estimate of survival at the corresponding time

# Kaplan-Meyer plot by ggsurvfit
survfit2(Surv(time, status) ~ 1, data = lung) |> 
  ggsurvfit() +
  labs(
    x = "Days",
    y = "Overall survival probability"
  ) +
add_confidence_interval() +
add_risktable()

# SURVIVAL: estimating-year survival--------------------------------------------

# Put 1 year as arbitrary endpoint (365.25 days)
summary(survfit(Surv(time, status) ~ 1, data = lung), times = 365.25)
# Output: survival 40,9%, incorrect naive estimation 47% if we add people died (censored) before 1 year of survival

#x-time survival probability estimates using the tbl_survfit()
survfit(Surv(time, status) ~ 1, data = lung) |> 
  tbl_survfit(
    times = 365.25,
    label_header = "**1-year survival (95% CI)**"
  )

# Estimating median survival time
survfit(Surv(time, status) ~ 1, data = lung)


# tables of median survival time estimates using the tbl_survfit() function from the {gtsummary} package:
survfit(Surv(time, status) ~ 1, data = lung) |> 
  tbl_survfit(
    probs = 0.5,
    label_header = "**Median survival (95% CI)**"
  )

# Comparing survival times between groups #-------------------------------------

# log-rank significance: comparing survival times between groups
survdiff(Surv(time, status) ~ sex, data = lung)


# The Cox regression mode #-----------------------------------------------------

# Quantify an effect size for a single variable, or multiple variables (coxph() from Surv())
coxph(Surv(time, status) ~ sex, data = lung)

# tables of results (tbl_regression() from {gtsummary})
# exponentiate set to TRUE to return the hazard ratio rather than the log hazard ratio
coxph(Surv(time, status) ~ sex, data = lung) |> 
  tbl_regression(exp = TRUE) 




# PART2: LANKDMARK ANALYSIS AND TIME DEPENDENT COVARIATES #---------------------

# BMT dataset #-----------------------------------------------------------------

# BMT dataset from {SemiCompRisks} package as an example dataset. 
# 137 bone marrow transplant patients. Variables of interest include:
# - T1 time (in days) to death or last follow-up;
# - delta1 death indicator; 1=Dead, 0=Alive;
# - TA time (in days) to acute graft-versus-host disease;
# - deltaA acute graft-versus-host disease indicator; 1-Developed acute graft-versus-host disease, 0-Never developed acute graft-versus-host disease


#install.packages("SemiCompRisks")
data(BMT, package = "SemiCompRisks")

head(BMT[, c("T1", "delta1", "TA", "deltaA")])

# Select a fixed time after baseline as your landmark time (90 days in this case).  
# Note: this should be done based on clinical information, prior to data inspection.

# Subset population for those followed at least until landmark time
lm_dat <- 
  BMT |> 
  filter(T1 >= 90) 

# 15 patients are excluded because died before 90 ds

# Calculate follow-up time from landmark and apply traditional methods
lm_dat <- 
  lm_dat |> 
  mutate(
    lm_T1 = T1 - 90
  )

survfit2(Surv(lm_T1, delta1) ~ deltaA, data = lm_dat) |> 
  ggsurvfit() +
  labs(
    x = "Days from 90-day landmark",
    y = "Overall survival probability"
  ) +
  add_risktable()

# coxph using subset to exclude died patients
coxph(
  Surv(T1, delta1) ~ deltaA, 
  subset = T1 >= 90, 
  data = BMT
) |> 
  tbl_regression(exp = TRUE)


# Alternative Time-dependent covariate approach --------------------------------

# The alternative more appropriate than landmark analysis when:
# - the value of a covariate is changing over time
# - there is not an obvious landmark time
# - use of a landmark would lead to many exclusions


# There was no ID variable in the BMT data, which is needed to create the special 
# dataset, so create an ID variable called my_id
BMT <- rowid_to_column(BMT, "my_id") #rowid_to_column in tibble package

# Use the tmerge function with the event and tdc function options to create the special dataset.
# - tmerge() creates a long dataset with multiple time intervals for the different covariate values for each patient
# - event() creates the new event indicator to go with the newly created time intervals
# - tdc() creates the time-dependent covariate indicator to go with the newly created time intervals

td_dat <- 
  tmerge(
    data1 = BMT |> select(my_id, T1, delta1), 
    data2 = BMT |> select(my_id, T1, delta1, TA, deltaA), 
    id = my_id, 
    death = event(T1, delta1),
    agvhd = tdc(TA)
  )


# Analyze this time-dependent covariate with coxph and an alteration to our use 
# of Surv to include arguments to both time and time2
coxph(
  Surv(time = tstart, time2 = tstop, event = death) ~ agvhd, 
  data = td_dat
) |> 
  tbl_regression(exp = TRUE)




# PART 3: COMPETING RISKS ------------------------------------------------------


# Melanoma dataset --------------------------------------------------------------

# From {MASS} package. Variables:
#- time survival time in days, possibly censored.
# - status 1 died from melanoma, 2 alive, 3 dead from other causes.
# - sex 1 = male, 0 = female.
# - age age in years.
# - year of operation.
# - thickness tumor thickness in mm.
# - ulcer 1 = presence, 0 = absence.

#install.packages("MASS")
data(Melanoma, package = "MASS")

# Recode status variable in a standard way (0=alive, 1=died from melanoma, 2=dead from other causes)
Melanoma <- 
  Melanoma |> 
  mutate(
    status = as.factor(case_when(status, `2` = 0, `1` = 1, `3` = 2))
  )

head(Melanoma)


# Cumulative incidence for competing risks
cuminc(Surv(time, status) ~ 1, data = Melanoma)

# By default ggcuminc() plots the first event type only. the following plot shows the cumulative incidence of death from melanoma
cuminc(Surv(time, status) ~ 1, data = Melanoma) |> 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()

# include both event types
cuminc(Surv(time, status) ~ 1, data = Melanoma) |> 
  ggcuminc(outcome = c("1", "2")) +
  ylim(c(0, 1)) + 
  labs(
    x = "Days"
  )

# Presence or absence of ulceration. Estimate the cumulative incidence at various times 
# by group and display that in a table using the tbl_cuminc() + Gray’s test to test for 
# a difference between groups over the entire follow-up period using the add_p() function.
cuminc(Surv(time, status) ~ ulcer, data = Melanoma) |> 
  tbl_cuminc(
    times = 1826.25, 
    label_header = "**{time/365.25}-year cuminc**") |> 
  add_p()


#plot of death due to melanoma, according to ulceration status
cuminc(Surv(time, status) ~ ulcer, data = Melanoma) |> 
  ggcuminc() + 
  labs(
    x = "Days"
  ) + 
  add_confidence_interval() +
  add_risktable()



# Competing risks regression

#crr() estimates the subdistribution hazards
crr(Surv(time, status) ~ sex + age, data = Melanoma)

#table tbl_regression(), with exp = TRUE to obtain the hazard ratio estimates
crr(Surv(time, status) ~ sex + age, data = Melanoma) |> 
  tbl_regression(exp = TRUE)

# OR in case need to censor death by other causes
coxph(
  Surv(time, ifelse(status == 1, 1, 0)) ~ sex + age, 
  data = Melanoma
) |> 
  tbl_regression(exp = TRUE)





# PART 4: ADVANCED TOPICS ------------------------------------------------------

# 1. Assessing the proportional hazards assumption
# 2. Making a smooth survival plot based on x-year survival according to a continuous covariate
# 3. Conditional survival


# 1. Assessing proportional hazards --------------------------------------------

# Assumption of the Cox prop. hazards regr. model is that the hazards are proportional at each point in time throughout follow-up (cox.zph() from {survival})

mv_fit <- coxph(Surv(time, status) ~ sex + age, data = lung)
cz <- cox.zph(mv_fit)
print(cz)

plot(cz)

# p-values >0.05 --> do not reject the null hypothesis, the proportional hazards assumption is satisfied for each individual covariate


# Smooth survival plot (continuous variable)
sm.options(
  list(
    xlab = "Age (years)",
    ylab = "Median time to death (years)")
)

sm.survival(
  x = lung$age,
  y = lung$time,
  status = lung$status,
  h = (1/6) * sd(lung$age) / nrow(lung)^(-1/4)
)


# Conditional survival: survival estimates among a group of patients who have already survived for some length of time
fit1 <- survfit(Surv(time, status) ~ 1, data = lung)

prob_times <- seq(365.25, 182.625 * 4, 182.625)

purrr::map_df(
  prob_times, 
  ~conditional_surv_est(
    basekm = fit1, 
    t1 = 182.625, 
    t2 = .x) 
) |> 
  mutate(months = round(prob_times / 30.4)) |> 
  select(months, everything()) |> 
  kable()


# visualize conditional survival data based on different lengths of time survived
gg_conditional_surv(
  basekm = fit1, 
  at = prob_times,
  main = "Conditional survival in lung data",
  xlab = "Days"
) +
  labs(color = "Conditional time")

