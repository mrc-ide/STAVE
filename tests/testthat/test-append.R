test_that("examples of append failing on invalid data", {
  
  # read in correctly formatted data from inst/extdata
  file_path <- system.file("extdata", "structure_tests/preload_01.xlsx", package = "STAVE")
  input_studies <- readxl::read_excel(file_path, sheet = "studies")
  input_surveys <- readxl::read_excel(file_path, sheet = "surveys")
  input_counts <- readxl::read_excel(file_path, sheet = "counts")
  
  # make object and load data
  z <- STAVE_object$new()
  z$append_data(studies_dataframe = input_studies,
                surveys_dataframe = input_surveys,
                counts_dataframe = input_counts) |>
    suppressMessages()
  
  # list of file names containing invalid data
  file_list <- sprintf("invalid_structure_0%s", 1:9)
  
  # check that all files in list fail
  for (i in seq_along(file_list)) {
    
    file_path <- system.file("extdata", "structure_tests", package = "STAVE")
    x_path <- sprintf("%s/%s.xlsx", file_path, file_list[i])
    x1 <- readxl::read_xlsx(x_path, sheet = "studies")
    x2 <- readxl::read_xlsx(x_path, sheet = "surveys")
    x3 <- readxl::read_xlsx(x_path, sheet = "counts")
    
    z$append_data(studies_dataframe = x1,
                  surveys_dataframe = x2,
                  counts_dataframe = x3) |>
      expect_error()
    
  }
  
})

test_that("examples of append working when data contains many missing fields", {
  
  # read in correctly formatted data from inst/extdata
  file_path <- system.file("extdata", "structure_tests/correct_01.xlsx", package = "STAVE")
  input_studies <- readxl::read_excel(file_path, sheet = "studies")
  input_surveys <- readxl::read_excel(file_path, sheet = "surveys")
  input_counts <- readxl::read_excel(file_path, sheet = "counts")
  
  # make object and load data
  z <- STAVE_object$new()
  expect_error(suppressMessages(z$append_data(studies_dataframe = input_studies,
                                              surveys_dataframe = input_surveys,
                                              counts_dataframe = input_counts)), NA)
  
})

test_that("further example of data that should pass", {
  
  # read in correctly formatted data from inst/extdata
  file_path <- system.file("extdata", "structure_tests/invalid_structure_10.xlsx", package = "STAVE")
  input_studies <- readxl::read_excel(file_path, sheet = "studies")
  input_surveys <- readxl::read_excel(file_path, sheet = "surveys")
  input_counts <- readxl::read_excel(file_path, sheet = "counts")
  
  # make object and load data
  z <- STAVE_object$new()
  expect_error(suppressMessages(z$append_data(studies_dataframe = input_studies,
                                              surveys_dataframe = input_surveys,
                                              counts_dataframe = input_counts)), NA)
  
})
