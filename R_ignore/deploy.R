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
# Useful commands:
# pkgdown::build_site() # build all pages of pkgdown website
# pkgdown::build_article('core_design_principles')  # build single vignette
# browseURL("docs/articles/core_design_principles.html")

# ------------------------------------------------------------------

data("example_input")
example_input


# create new object
s <- STAVE_object$new()
s

private <- s$.__enclos_env__$private

# append pre-loaded data
s$append_data(studies_dataframe = example_input$studies,
              surveys_dataframe = example_input$surveys,
              counts_dataframe = example_input$counts)
s

# calculate prevalence
p <- s$get_prevalence(target_variant = "k13:469:Y",
                      keep_ambiguous = TRUE,
                      prev_from_min = TRUE,
                      return_full = TRUE)
p

# get list of variant loci
s$get_variants()

# drop a study
d <- s$get_studies()$study_id
s$drop_study(drop_study_id = d[2])
s$drop_survey(drop_survey_id = c("Asua_2019_Agago_2017"))
s
