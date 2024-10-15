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

# get list of variant loci
s$get_variants()

# drop a study
d <- s$get_studies()$study_ID
s$drop_study(drop_study_ID = d[1:5])
s
