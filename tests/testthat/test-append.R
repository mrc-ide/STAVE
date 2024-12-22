test_that("test that fails when appending data that clashes with existing preloaded data", {
  
  # list of file names containing invalid data
  file_list <- sprintf("invalid_append_%02d", 1:1)
  
  # check that all files in list fail
  for (i in seq_along(file_list)) {
    
    folder_path <- system.file("extdata", "append_tests", package = "STAVE")
    preload_path <- sprintf("%s/%s_preload.xlsx", folder_path, file_list[i])
    file_path <- sprintf("%s/%s.xlsx", folder_path, file_list[i])
    
    # preload correct data
    preload_studies <- readxl::read_xlsx(preload_path, sheet = "studies")
    preload_surveys <- readxl::read_xlsx(preload_path, sheet = "surveys")
    preload_counts <- readxl::read_xlsx(preload_path, sheet = "counts")
    
    z <- STAVE_object$new()
    z$append_data(studies_dataframe = preload_studies,
                  surveys_dataframe = preload_surveys,
                  counts_dataframe = preload_counts) |>
      suppressMessages()
    
    # try to append bad data
    input_studies <- readxl::read_xlsx(file_path, sheet = "studies")
    input_surveys <- readxl::read_xlsx(file_path, sheet = "surveys")
    input_counts <- readxl::read_xlsx(file_path, sheet = "counts")
    
    z$append_data(studies_dataframe = input_studies,
                  surveys_dataframe = input_surveys,
                  counts_dataframe = input_counts) |>
      suppressMessages() |>
      expect_error()
  }
  
})
