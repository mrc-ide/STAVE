test_that("get_prevalence() checks calculations return the expected values", {
  
  # read in correctly formatted data from inst/extdata
  file_path <- system.file("extdata", "prevalence_tests/correct_01.xlsx", package = "STAVE")
  input_studies <- readxl::read_excel(file_path, sheet = "studies")
  input_surveys <- readxl::read_excel(file_path, sheet = "surveys")
  input_counts <- readxl::read_excel(file_path, sheet = "counts")
  input_prevalence <- readxl::read_excel(file_path, sheet = "prevalence")
  input_prevalence2 <- readxl::read_excel(file_path, sheet = "prevalence2")
  
  # make object and load data
  z <- STAVE_object$new()
  z$append_data(studies_dataframe = input_studies,
                surveys_dataframe = input_surveys,
                counts_dataframe = input_counts) |>
    suppressMessages()
  
  # test all prevalences in the input_prevalence table
  for (i in 1:nrow(input_prevalence)) {
    
    # calculate prevalence
    p <- z$get_prevalence(target_variant = input_prevalence$target_variant[i],
                          keep_ambiguous = input_prevalence$keep_ambiguous[i],
                          prev_from_min = input_prevalence$prev_from_min[i]) |>
      suppressWarnings()
    
    expect_equal(p$numerator, input_prevalence$numerator[i])
    expect_equal(p$denominator, input_prevalence$denominator[i])
  }
  
  # test all prevalences in the input_prevalence2 table
  for (i in 1:nrow(input_prevalence2)) {
    
    # calculate prevalence
    p <- z$get_prevalence(target_variant = input_prevalence$target_variant[i],
                          keep_ambiguous = input_prevalence$keep_ambiguous[i],
                          prev_from_min = input_prevalence$prev_from_min[i]) |>
      suppressWarnings()
    
    expect_equal(p$numerator, input_prevalence$numerator[i])
    expect_equal(p$denominator, input_prevalence$denominator[i])
  }
  
})
