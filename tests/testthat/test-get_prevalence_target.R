test_that("get_prevalence() target variant formatted correctly", {
  
  # read in correctly formatted data from inst/extdata
  file_path <- system.file("extdata", "prevalence_tests/correct_01.xlsx", package = "STAVE")
  input_studies <- readxl::read_excel(file_path, sheet = "studies")
  input_surveys <- readxl::read_excel(file_path, sheet = "surveys") |>
    mutate(collection_start = as.Date(collection_start),
           collection_end = as.Date(collection_end),
           collection_day = as.Date(collection_day))
  input_counts <- readxl::read_excel(file_path, sheet = "counts")
  
  # make object and load data
  z <- STAVE_object$new()
  z$append_data(studies_dataframe = input_studies,
                surveys_dataframe = input_surveys,
                counts_dataframe = input_counts) |>
    suppressMessages()
  
  # get prevalence with badly formatted targets
  expect_error(suppressMessages(z$get_prevalence(target_variant = "crt")))
  expect_error(z$get_prevalence(target_variant = "crt:1:A/C"))
  expect_error(z$get_prevalence(target_variant = "crt:1:A|C"))
  
})
