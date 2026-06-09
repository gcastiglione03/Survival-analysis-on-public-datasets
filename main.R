# ==============================================================================
# main.R
# Run all analysis scripts in order
# ==============================================================================

library(here)

source(here("R", "01_kaplan_meier.R"))
source(here("R", "02_landmark_time_dependent.R"))
source(here("R", "03_competing_risks.R"))
source(here("R", "04_advanced_topics.R"))