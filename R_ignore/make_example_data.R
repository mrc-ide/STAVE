# make_example_data.R
#
# Author: Bob Verity
# Date: 2025-11-29
#
# Inputs: R_ignore/example_data.xlsx
#
# Outputs: data/example_inputs.rda
#
# Purpose:
# Reads in example data, cleans and saves as rda.
#
# ------------------------------------------------------------------

library(tidyverse)
library(here)
library(readxl)

# NOTE, advised to clear session before running, as existing data loaded from
# package seems to be interfering in some cases
#rm(list = ls())

# read data from file
df_studies <- readxl::read_excel(here("R_ignore", "example_data.xlsx"), sheet = "studies")
df_surveys <- readxl::read_excel(here("R_ignore", "example_data.xlsx"), sheet = "surveys")
df_counts <- readxl::read_excel(here("R_ignore", "example_data.xlsx"), sheet = "counts")

# fix some formats
df_studies <- df_studies |>
  mutate(description = as.character(description))

df_surveys <- df_surveys |>
  mutate(collection_start = as.Date(collection_start),
         collection_end = as.Date(collection_end),
         collection_day = as.Date(collection_day),
         location_notes = as.character(location_notes),
         time_notes = as.character(time_notes))

# make combined object and save to file
example_input <- list(studies = df_studies,
                      surveys = df_surveys,
                      counts = df_counts)

save(example_input, file = here("data", "example_input.rda"))
