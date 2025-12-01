test_that("tests of failing to load bad data due to structural problems", {
  
  # list of file names containing invalid data
  file_list <- sprintf("invalid_structure_%02d", 1:10)
  
  # check that all files in list fail
  for (i in seq_along(file_list)) {
    
    folder_path <- system.file("extdata", "structure_tests", package = "STAVE")
    file_path <- sprintf("%s/%s.xlsx", folder_path, file_list[i])
    input_studies <- readxl::read_xlsx(file_path, sheet = "studies")
    input_surveys <- readxl::read_xlsx(file_path, sheet = "surveys") |>
      mutate(collection_start = as.Date(collection_start),
             collection_end = as.Date(collection_end),
             collection_day = as.Date(collection_day))
    input_counts <- readxl::read_xlsx(file_path, sheet = "counts")
    
    z <- STAVE_object$new()
    z$append_data(studies_dataframe = input_studies,
                  surveys_dataframe = input_surveys,
                  counts_dataframe = input_counts) |>
      suppressMessages() |>
      expect_error()
  }
})

test_that("tests of successfully loading data even when there is complex structure", {
  
  # list of file names containing invalid data
  file_list <- sprintf("correct_structure_%02d", 1)
  
  # check that all files in list pass
  for (i in seq_along(file_list)) {
    
    folder_path <- system.file("extdata", "structure_tests", package = "STAVE")
    file_path <- sprintf("%s/%s.xlsx", folder_path, file_list[i])
    input_studies <- readxl::read_xlsx(file_path, sheet = "studies")
    input_surveys <- readxl::read_xlsx(file_path, sheet = "surveys") |>
      mutate(collection_start = as.Date(collection_start),
             collection_end = as.Date(collection_end),
             collection_day = as.Date(collection_day))
    input_counts <- readxl::read_xlsx(file_path, sheet = "counts")
    
    z <- STAVE_object$new()
    z$append_data(studies_dataframe = input_studies,
                  surveys_dataframe = input_surveys,
                  counts_dataframe = input_counts) |>
      suppressMessages() |>
      expect_error(NA)
  }
})

