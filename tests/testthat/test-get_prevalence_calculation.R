test_that("get_prevalence() checks calculations return the expected values", {
  
  # read in correctly formatted data from inst/extdata
  file_path <- system.file("extdata", "prevalence_tests/correct_01.xlsx", package = "STAVE")
  input_studies <- readxl::read_excel(file_path, sheet = "studies")
  input_surveys <- readxl::read_excel(file_path, sheet = "surveys")
  input_counts <- readxl::read_excel(file_path, sheet = "counts")
  
  # make object and load data
  z <- STAVE_object$new()
  z$append_data(studies_dataframe = input_studies,
                surveys_dataframe = input_surveys,
                counts_dataframe = input_counts) |>
    suppressMessages()
  
  # calculate prevalence
  p <- z$get_prevalence("crt:72:C")
  expect_equal(p$numerator, 4)
  expect_equal(p$denominator, 30)
  
  p <- z$get_prevalence("crt:72_73:A_A", keep_ambiguous = TRUE, prev_from_min = TRUE)
  expect_equal(p$numerator_min, 1)
  expect_equal(p$numerator_max, 2)
  expect_equal(p$denominator, 20)
  
  p <- z$get_prevalence("k13:580:Y")
  expect_equal(p$numerator, 0)
  expect_equal(p$denominator, 0)
  
})
