
#' @title STAVE data object (R6 class)
#'
#' @description 
#' TODO
#' 
#' @importFrom R6 R6Class
#' @export

STAVE_object <- R6::R6Class(
  "STAVE_object",
  
  # Public fields and methods
  public = list(
    
    # -----------------------------------
    #' @description
    #' Custom print method to control console output
    print = function() {
      cat("TODO - custom print")
    },
    
    # -----------------------------------
    #' @description
    #' Append new data
    #' @param studies_dataframe TODO
    #' @param surveys_dataframe TODO
    #' @param counts_dataframe TODO
    append_data = function(studies_dataframe, surveys_dataframe, counts_dataframe) {
      private$counts <- counts_dataframe
    },
    
    # -----------------------------------
    #' @description
    #' Calculate prevalence
    #' @param variant TODO
    get_prevalence = function(variant) {
      cat("TODO - calculate prevalence")
    },
    
    # -----------------------------------
    #' @description
    #' Return a table of all variants present in the data object.
    #' @param report_haplo (Boolean) if TRUE then list all haplotypes.
    #'   Otherwise, list in locus-by-locus format. Defaults to FALSE.
    get_variants = function(report_haplo = FALSE) {
      cat("TODO - return table of variants")
    },
    
    # -----------------------------------
    #' @description
    #' Export full data object as an Excel spreadsheet (.xlsx format).
    #' @param file_path the full path and name of the output file.
    export_xlsx = function(file_path) {
      cat("TODO - export to xlsx")
    }
    
  ),
  
  
  # Private fields and methods
  private = list(
    
    # main tables
    studies = NULL,
    surveys = NULL,
    counts = NULL,
    
    # -----------------------------------
    # test
    show_secret = function() {
      cat("Secret Message:", private$counts, "\n")
    }
  )
)
