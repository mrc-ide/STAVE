# deploy.R
#
# Author: Bob Verity
# Date: 2024-10-08
#
# Inputs: (none)
#
# Outputs: (none)
#
# Purpose:
# Example deploy script for STAVE package.
#
# ------------------------------------------------------------------

# create new object
s <- STAVE_object$new()
s

# append pre-loaded data
s$append_data(studies_dataframe = example_input$studies,
              surveys_dataframe = example_input$surveys,
              counts_dataframe = example_input$counts)
s

# calculate prevalence
p <- s$get_prevalence("mdr1:184:F")
p

s$get_variants()


v <- s$get_counts()$variant_string

v_split <- strsplit(v, split = ";") |>
  unlist()


v_split
