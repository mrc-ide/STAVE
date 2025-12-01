
#------------------------------------------------
#' @title STAVE data object (R6 class)
#'
#' @description 
#' The main class that stores the data and is responsible for all data input,
#' output, and processing functions. Most of the functionality of the STAVE
#' package is through this class in the form of member functions.
#' 
#' @details
#' The raw data are stored as private variables within this object, meaning they
#' cannot (or should not) be edited directly. Rather, tables can be extracted
#' using \code{get_counts()} and similarly for other tables. The three tables
#' are:
#'   \enumerate{
#'     \item studies: Information on where the data came from, for example a url
#'     and author names. Each study is indexed with a unique study_id.
#'     \item surveys: Information on the surveys represented within a study. A
#'     survey is defined here as a discrete instance of data collection, which
#'     includes information on geography (latitude and longitude) and collection
#'     time. Surveys are given survey_ids and are linked to a particular study
#'     through the study_id.
#'     \item counts: The actual genetic information, which is linked to a
#'     particular study and survey through the study_id and survey_id. Genetic
#'     variants are encoded in character strings that must follow a specified
#'     format, and the number of times this variant was observed among the total
#'     sample is stored in columns.
#'   }
#' This combination of linked tables allows efficient and flexible encoding of
#' variants, while avoiding unnecessary duplication of information.
#' 
#' @importFrom R6 R6Class
#' @import dplyr variantstring
#' @importFrom tidyr replace_na
#' @export

STAVE_object <- R6::R6Class(
  "STAVE_object",
  
  # Public fields and methods
  public = list(
    
    # -----------------------------------
    #' @description
    #' Custom print method to control console output
    print = function() {
      if (is.null(private$studies)) {
        cat("No data loaded")
      } else {
        if (nrow(private$studies) == 0) {
          cat("No data loaded")
        } else {
          cat(sprintf("Studies: %s", nrow(private$studies)))
          n_surveys <- length(unique(private$surveys$survey_id))
          cat(sprintf("\nSurveys: %s", n_surveys))
        }
      }
    },
    
    # -----------------------------------
    #' @description
    #' Extract the studies data.frames stored within the object
    get_studies = function() {
      return(private$studies)
    },
    
    # -----------------------------------
    #' @description
    #' Extract the surveys data.frames stored within the object
    get_surveys = function() {
      return(private$surveys)
    },
    
    # -----------------------------------
    #' @description
    #' Extract the counts data.frames stored within the object
    get_counts = function() {
      return(private$counts)
    },
    
    # -----------------------------------
    #' @description
    #' Extract the version number of the STAVE object. This is important as
    #' member functions of a STAVE object are directly linked to the object
    #' itself, and will not be updated by updating the version of the package in
    #' your environment. To update a STAVE object to a new package version, you
    #' should first extract the data and then load into a new STAVE object
    #' created with the most recent version.
    get_version = function() {
      return(private$version)
    },
    
    # -----------------------------------
    #' @description
    #' Append new data
    #' @param studies_dataframe a data frame containing information at the study
    #'   level. This data frame must have the following columns: study_id,
    #'   study_label, description, access_level, contributors, reference, reference_year,
    #'   notes. Compulsory fields (no missing values) are: study_id and reference.
    #' @param surveys_dataframe a data frame containing information at the
    #'   survey level. This data.frame must have the following columns:
    #'   study_id, survey_id, country_name, site_name, latitude, longitude, location_method,
    #'   location_notes, collection_start, collection_end, collection_day, time_method,
    #'   time_notes. Compulsory fields (no missing values) are: study_id, survey_id,
    #'   latitude, longitude, collection_day.
    #' @param counts_dataframe a data.frame of genetic information. Must contain
    #'   the following columns: study_id, survey_id, variant_string, variant_num,
    #'   total_num. All fields are compulsory (no missing values).
    append_data = function(studies_dataframe, surveys_dataframe, counts_dataframe) {
      
      # basic checks on format of input data.frames. These checks ensure that
      # each independent row of the input data are correctly formatted
      private$check_format_studies(studies_dataframe)
      private$check_format_surveys(surveys_dataframe)
      private$check_format_counts(counts_dataframe)
      
      # structural checks on data. These are more focused on the relationships
      # between data, for example ensuring no duplicate study_ids, or ensuring
      # that the same variant codons are not encoded multiple times within a
      # single survey
      private$check_structure(studies_dataframe, surveys_dataframe, counts_dataframe)
      
      # append to existing tables
      private$studies <- rbind(private$studies, studies_dataframe)
      private$surveys <- rbind(private$surveys, surveys_dataframe)
      private$counts <- rbind(private$counts, counts_dataframe)
      
      message("data correctly appended")
    },
    
    # -----------------------------------
    #' @description
    #' Calculate prevalence
    #' @param target_variant the variant on which to calculate prevalence, for
    #'   example crt:72:C. There can be no heterozygous calls within this
    #'   string.
    #' @param keep_ambiguous there may be variants in the data for which the
    #'   target_variant could be in the sample, but this cannot be proven
    #'   conclusively. For example, the sequence A_A_A may be a match to the
    #'   sequence A/C_A/C_A, or it may not. These are unphased genotypes so we
    #'   cannot be sure. If \code{keep_ambiguous = TRUE} then both a min and a
    #'   max numerator are reported that either include all ambiguous calls as
    #'   matches (max) or exclude them as mismatches (min). If \code{FALSE} (the
    #'   default) then ambiguous calls are skipped over, which may downwardly
    #'   bias prevalence calculation.
    #' @param prev_from_min the output object includes a point estimate of the
    #'   prevalence along with exact binomial confidence intervals. In the case
    #'   of ambiguous calls, these must be calculated from one of
    #'   \code{numerator_min} or \code{numerator_max}. This argument sets which
    #'   one of these values is used in the calculation. Defaults to
    #'   \code{TRUE}, which risks underestimating prevalence (whereas the
    #'   alternative risks overestimating prevalence).
    #' @param return_full if \code{TRUE} (the default) returns the entire loaded
    #'   dataset, with prevalence equal to \code{NA} if there is no denominator.
    #'   If \code{FALSE} only returns entries for which there is a non-zero
    #'   denominator.
    #'   
    #'   @import dplyr
    #'   @importFrom tidyr replace_na
    get_prevalence = function(target_variant, keep_ambiguous = FALSE, prev_from_min = TRUE, return_full = TRUE) {
      
      # basic input checks
      assert_single_string(target_variant)
      assert_single_logical(keep_ambiguous)
      assert_single_logical(prev_from_min)
      variantstring::check_variant_string(target_variant)
      
      # extract target position
      target_position <- variantstring::position_from_variant_string(target_variant)
      
      # compare variant strings to get matches
      df_match <- variantstring::compare_variant_string(target_string = target_variant,
                                                        comparison_strings = private$counts$variant_string)
      
      df_count_match <- private$counts |>
        mutate(match = df_match$match,
               ambiguous = df_match$ambiguous)
      
      # inform user if there are ambiguous matches to consider
      any_ambiguous <- any(df_count_match$ambiguous)
      if ((keep_ambiguous == FALSE) & any_ambiguous) {
        warning(wrap_message(paste("There are some ambiguous matches to the target sequence due to mixed",
                                   "infections. Consider running with keep_ambiguous = TRUE to return both",
                                   "upper and lower bounds on the numerator. Currently ignoring all ambiguous",
                                   "matches, which may lead to prevalence being downwardly biased.",
                                   collapse = " ")))
      }
      
      # look for match based on position, used in denominator calculation
      df_count_match <- df_count_match |>
        mutate(position_string = variantstring::position_from_variant_string(variant_string),
               match_pos = variantstring::compare_position_string(
                 target_string = target_position,
                 comparison_strings = variant_string))
      
      # calculate numerator and denominator of prevalence calculation
      df_prev <- df_count_match |>
        filter(match_pos) |>
        group_by(study_id, survey_id, position_string) |>
        summarise(numerator_min = sum(variant_num * match * (1 - ambiguous)),
                  numerator_max = sum(variant_num * match),
                  denominator = total_num[1],
                  .groups = "drop") |>
        ungroup() |>
        select(-position_string) |>
        group_by(study_id, survey_id) |>
        summarise(numerator_min = sum(numerator_min),
                  numerator_max = sum(numerator_max),
                  denominator = sum(denominator),
                  .groups = "drop")
      
      # calculate 95% exact binomial CI
      df_prev <- df_prev |>
        mutate(numerator = if(prev_from_min) numerator_min else numerator_max,
               prevalence = 100 * numerator / denominator,
               prevalence_lower = 100 * qbeta(0.025, shape1 = numerator, shape2 = denominator - numerator + 1),
               prevalence_upper = 100 * qbeta(0.975, shape1 = numerator + 1, shape2 = denominator - numerator)) |>
        relocate(numerator, .before = numerator_min)
      
      # tidy up names
      if (!keep_ambiguous) {
        df_prev <- df_prev |>
          select(-numerator_min, -numerator_max)
      }
      
      # merge back with all study and survey info
      ret <- private$studies |>
        left_join(private$surveys, by = join_by(study_id)) |>
        left_join(df_prev, by = join_by(study_id, survey_id)) |>
        mutate(numerator = tidyr::replace_na(numerator, 0),
               denominator = tidyr::replace_na(denominator, 0))
      
      if (keep_ambiguous) {
        ret <- ret |>
          mutate(numerator_min = tidyr::replace_na(numerator_min, 0),
                 numerator_max = tidyr::replace_na(numerator_max, 0))
      }
      
      # optional filter out 0s
      if (!return_full) {
        ret <- ret |>
          filter(denominator > 0)
      }
      
      return(ret)
    },
    
    # -----------------------------------
    #' @description
    #' Return a vector of all variants present in the data object.
    #' @param report_haplo (Boolean) if TRUE then list all haplotypes.
    #'   Otherwise, list in locus-by-locus format. Defaults to FALSE.
    get_variants = function(report_haplo = FALSE) {
      
      v <- self$get_counts()$variant_string
      if (is.null(v)) {
        stop("No data loaded")
      }
      
      # simple option to return all unique haplos
      if (report_haplo) {
        ret <- sort(unique(v))
        return(ret)
      }
      
      # split into per-locus variants
      ret <- variantstring::extract_single_locus_variants(v) |>
        unlist() |>
        unique() |>
        sort()
      
      return(ret)
    },
    
    # -----------------------------------
    #' @description
    #' Drop one or more study_ids from the data. This will drop from all
    #' internally stored data objects, including the corresponding surveys and
    #' counts data.
    #' @param drop_study_id a vector of study_ids to drop from all data objects.
    drop_study = function(drop_study_id) {
      
      # check that study is a valid ID
      df_studies <- self$get_studies()
      if (!all(drop_study_id %in% df_studies$study_id)) {
        stop("all values in drop_study_id must be present in the loaded data")
      }
      
      # get corresponding survey IDs
      drop_survey_id <- self$get_surveys() |>
        select(study_id, survey_id) |>
        filter(study_id %in% drop_study_id)
      
      # drop IDs and save back into private objects
      private$studies <- self$get_studies() |>
        filter(!(study_id %in% drop_study_id))
      
      private$surveys <- self$get_surveys() |>
        filter(!(study_id %in% drop_study_id))
      
      private$counts <- self$get_counts() |>
        anti_join(drop_survey_id, by = join_by(study_id, survey_id))
      
      n_study_drop <- length(drop_study_id)
      n_survey_drop <- length(unique(drop_survey_id$survey_id))
      message(sprintf("drop %s %s, %s %s",
                      n_study_drop,
                      ifelse(n_study_drop == 1, "study", "studies"),
                      n_survey_drop,
                      ifelse(n_survey_drop == 1, "survey", "surveys")))
    },
    
    # -----------------------------------
    #' @description
    #' Drop one or more survey_ids from the data. This will drop from all
    #' internally stored data objects, including the corresponding counts data,
    #' and the study if this is the only survey in that study.
    #' @param drop_survey_id a vector of survey_ids to drop from all data objects.
    drop_survey = function(drop_survey_id) {
      
      # check that survey is a valid ID
      df_surveys <- self$get_surveys()
      if (!all(drop_survey_id %in% df_surveys$survey_id)) {
        stop("all values in drop_survey_id must be present in the loaded data")
      }
      
      # get number of surveys that need to be dropped from each study
      df_study_drop <- df_surveys |>
        filter(survey_id %in% drop_survey_id) |>
        group_by(study_id) |>
        summarise(n_drop = n())
      
      # get total surveys per study and merge to see what proportion are being dropped
      df_study_total <- df_surveys |>
        filter(study_id %in% df_study_drop$study_id) |>
        group_by(study_id) |>
        summarise(n_total = n()) |>
        ungroup() |>
        left_join(df_study_drop, by = join_by(study_id))
      
      # get study IDs to drop
      drop_study_id <- df_study_total |>
        filter(n_total == n_drop) |>
        pull(study_id)
      
      # drop IDs and save back into private objects
      if (length(drop_study_id) > 0) {
        private$studies <- self$get_studies() |>
          filter(!(study_id %in% drop_study_id))
      }
      
      private$surveys <- df_surveys |>
        filter(!(survey_id %in% drop_survey_id))
      
      private$counts <- self$get_counts() |>
        filter(!(survey_id %in% drop_survey_id))
      
      n_study_drop <- length(drop_study_id)
      n_survey_drop <- length(drop_survey_id)
      message(sprintf("drop %s %s, %s %s",
                      n_study_drop,
                      ifelse(n_study_drop == 1, "study", "studies"),
                      n_survey_drop,
                      ifelse(n_survey_drop == 1, "survey", "surveys")))
    }
  ),
  
  # Private fields and methods
  private = list(
    
    # main tables
    studies = NULL,
    surveys = NULL,
    counts = NULL,
    version = packageVersion("STAVE"),
    
    # -----------------------------------
    # check basic formatting (variables types etc) of studies data.frame
    check_format_studies = function(studies_dataframe) {
      
      # CHECKS
      # - must be data.frame with exactly the right number of columns and correctly named
      # - variables have the correct type and range. Allow for missing data in some cases
      # - IDs are valid identifiers
      # - warning if any studies have study_type of 'private'
      
      # basic structure
      assert_dataframe(studies_dataframe)
      assert_ncol(studies_dataframe, 7)
      assert_eq(names(studies_dataframe), c("study_id", "study_label", "description", "access_level",
                                            "contributors", "reference", "reference_year"))
      
      # study_id (compulsory)
      assert_non_NA(studies_dataframe$study_id)
      assert_valid_string(studies_dataframe$study_id, message_name = "study_id in studies_dataframe")
      
      # study_label
      w <- which(!is.na(studies_dataframe$study_label))
      if (any(w)) {
        assert_string(studies_dataframe$study_label[w])
      }
      
      # study_label
      w <- which(!is.na(studies_dataframe$description))
      if (any(w)) {
        assert_string(studies_dataframe$description[w])
      }
      
      # access_level
      assert_non_NA(studies_dataframe$access_level)
      assert_in(studies_dataframe$access_level, c("public", "restricted", "private"))
      if (any(studies_dataframe$access_level %in% c("restricted", "private"))) {
        warning(paste0("Some appended data are labelled as restricted or private.",
                       " Be sure to remove these before making this STAVE object public."))
      }
      
      # contributors
      w <- which(!is.na(studies_dataframe$contributors))
      if (any(w)) {
        assert_string(studies_dataframe$contributors[w])
      }
      
      # reference (compulsory)
      assert_non_NA(studies_dataframe$reference)
      assert_string(studies_dataframe$reference)
      
      # reference_year
      w <- which(!is.na(studies_dataframe$reference_year))
      if (any(w)) {
        assert_int(studies_dataframe$reference_year[w])
      }
      
      invisible(TRUE)
    },
    
    # -----------------------------------
    # check basic formatting (variables types etc) of surveys data.frame
    check_format_surveys = function(surveys_dataframe) {
      
      # CHECKS
      # - must be data.frame with exactly the right number of columns and correctly named
      # - variables have the correct type and range. Allow for missing data in some cases
      # - IDs are valid identifiers
      # - collection days are valid date strings
      
      # basic structure
      assert_dataframe(surveys_dataframe)
      assert_ncol(surveys_dataframe, 13)
      assert_eq(names(surveys_dataframe), c("study_id", "survey_id", "country_name", "site_name", "latitude",
                                            "longitude", "location_method", "location_notes", "collection_start",
                                            "collection_end", "collection_day", "time_method", "time_notes"))
      
      # study_key (compulsory)
      assert_non_NA(surveys_dataframe$study_id)
      assert_valid_string(surveys_dataframe$study_id, message_name = "study_id in surveys_dataframe")
      
      # survey_id (compulsory)
      assert_non_NA(surveys_dataframe$survey_id)
      assert_valid_string(surveys_dataframe$survey_id, message_name = "survey_id in surveys_dataframe")
      
      # country_name
      w <- which(!is.na(surveys_dataframe$country_name))
      if (any(w)) {
        assert_string(surveys_dataframe$country_name[w])
      }
      
      # site name
      w <- which(!is.na(surveys_dataframe$site_name))
      if (any(w)) {
        assert_string(surveys_dataframe$site_name[w])
      }
      
      # latitude (compulsory)
      assert_non_NA(surveys_dataframe$latitude)
      assert_bounded(surveys_dataframe$latitude, left = -90, right = 90)
      
      # longitude (compulsory)
      assert_non_NA(surveys_dataframe$longitude)
      assert_bounded(surveys_dataframe$longitude, left = -180, right = 180)
      
      # location_method
      w <- which(!is.na(surveys_dataframe$location_method))
      if (any(w)) {
        assert_string(surveys_dataframe$location_method[w])
      }
      
      # location_notes
      w <- which(!is.na(surveys_dataframe$location_notes))
      if (any(w)) {
        assert_string(surveys_dataframe$location_notes[w])
      }
      
      # collection_start
      w <- which(!is.na(surveys_dataframe$collection_start))
      if (any(w)) {
        assert_class(surveys_dataframe$collection_start[w], "Date")
      }
      
      # collection_end
      w <- which(!is.na(surveys_dataframe$collection_end))
      if (any(w)) {
        assert_class(surveys_dataframe$collection_end[w], "Date")
      }
      
      # collection_day (compulsory)
      assert_non_NA(surveys_dataframe$collection_day)
      assert_class(surveys_dataframe$collection_day[w], "Date")
      
      # time_method
      w <- which(!is.na(surveys_dataframe$time_method))
      if (any(w)) {
        assert_string(surveys_dataframe$time_method[w])
      }
      
      # time_notes
      w <- which(!is.na(surveys_dataframe$time_notes))
      if (any(w)) {
        assert_string(surveys_dataframe$time_notes[w])
      }
      
      invisible(TRUE)
    },
    
    # -----------------------------------
    # check basic formatting (variables types etc) of counts data.frame
    check_format_counts = function(counts_dataframe) {
      
      # CHECKS
      # - must be data.frame with exactly the right number of columns and correctly named
      # - variables have the correct type and range. Allow for missing data in some cases
      # - IDs are valid identifiers
      # - check valid variant string
      
      # basic structure
      assert_dataframe(counts_dataframe)
      assert_ncol(counts_dataframe, 5)
      assert_eq(names(counts_dataframe), c("study_id", "survey_id", "variant_string",
                                           "variant_num", "total_num"))
      
      # study_id (compulsory)
      assert_non_NA(counts_dataframe$study_id)
      assert_valid_string(counts_dataframe$study_id, message_name = "study_id in counts_dataframe")
      
      # survey_id (compulsory)
      assert_non_NA(counts_dataframe$survey_id)
      assert_valid_string(counts_dataframe$survey_id, message_name = "survey_id in counts_dataframe")
      
      # variant_string (compulsory)
      assert_non_NA(counts_dataframe$variant_string)
      variantstring::check_variant_string(counts_dataframe$variant_string)
      
      # variant_num (compulsory)
      assert_non_NA(counts_dataframe$variant_num)
      assert_pos_int(counts_dataframe$variant_num, zero_allowed = TRUE)
      
      # total_num (compulsory)
      assert_non_NA(counts_dataframe$total_num)
      assert_pos_int(counts_dataframe$total_num, zero_allowed = FALSE)
      
      invisible(TRUE)
    },
    
    # -----------------------------------
    # check that new inputs can be appended without breaking overall structure
    check_structure = function(studies_dataframe, surveys_dataframe, counts_dataframe) {
      
      # CHECKS
      # studies_dataframe:
      # - study_ids are unique
      # - study_ids are not already present in pre-loaded study_ids
      # - study_ids are referenced at least once in surveys table
      
      # surveys_dataframe:
      # - survey_ids are unique (completely unique, not just unique within studies)
      # - study_ids can be found in studies table
      # - survey_ids are referenced at least once in counts table
      
      # counts_dataframe:
      # - combination of study_id and survey_id can be found in surveys table
      # - variant_strings are unique within each study-survey (even if listed in different order)
      # - total_num is consistent between variants that correspond to the same positions
      # - after splitting by study-survey and also by gene-locus combination (i.e. the first two characteristics
      #   of each encoded gene, but not by the observed amino acids):
      #       - variant_num cannot sum to more than total_num
      
      # studies_dataframe:
      assert_noduplicates(studies_dataframe$study_id)
      if (any(studies_dataframe$study_id %in% private$studies$study_id)) {
        stop(wrap_message(paste0("studies_dataframe cannot contain any study_ids that already exist in the loaded",
                                 " studies table. If you want to add to an existing study then you must first",
                                 " drop the currently loaded version. See ?drop_study() for how to do this")))
      }
      if (!all(studies_dataframe$study_id %in% surveys_dataframe$study_id)) {
        stop(wrap_message(paste0("every study_id in the studies_dataframe must be referenced at least once in",
                                 " the surveys_dataframe")))
      }
      
      # surveys_dataframe:
      assert_noduplicates(surveys_dataframe$survey_id)
      if (!all(surveys_dataframe$study_id %in% studies_dataframe$study_id)) {
        stop(wrap_message("every study_id in the surveys_dataframe must be present in the studies_dataframe"))
      }
      if (!all(surveys_dataframe$survey_id %in% counts_dataframe$survey_id)) {
        stop(wrap_message(paste0("every survey_id in the surveys_dataframe must be referenced at least once in",
                                " the counts_dataframe")))
      }
      
      # counts_dataframe:
      counts_combined_id <- paste(counts_dataframe$study_id, counts_dataframe$survey_id, sep = ".")
      surveys_combined_id <- paste(surveys_dataframe$study_id, surveys_dataframe$survey_id, sep = ".")
      if (!all(counts_combined_id %in% surveys_combined_id)) {
        stop(wrap_message(paste0("every study_id-survey_id combination in the counts_dataframe must be present",
                                 " in the surveys_dataframe")))
      }
      
      # ensure no variant is present more than once (should not be after sorting alphabetically)
      variant_duplicated <- counts_dataframe |>
        mutate(variant_string_order = variantstring::order_variant_string(variant_string)) |>
        group_by(study_id, survey_id, variant_string_order) |>
        summarise(variant_duplicated = any(duplicated(variant_string_order)), .groups = "drop") |>
        pull(variant_duplicated)
      
      if (any(variant_duplicated)) {
        stop(wrap_message(paste0("the exact same variant cannot be present more than once in the same survey.",
                                 " This includes the same variant with the genes listed in different order")))
      }
      
      # ensure total_num is consistent between variants that correspond to the same positions
      distinct_total_num <- counts_dataframe |>
        mutate(position_string = variantstring::position_from_variant_string(variant_string)) |>
        group_by(study_id, survey_id, position_string) |>
        summarise(distinct_total_num = length(unique(total_num)), .groups = "drop") |>
        pull(distinct_total_num)
      
      if (any(distinct_total_num > 1)) {
        stop(wrap_message(paste0("if there are multiple rows in a survey that correspond to the same gene-locus",
                                 " combination then they must have the same total_num. This includes variants",
                                 " with heterozygous calls")))
      }
      
      # sum variants and determine if any where counts exceed the total_num
      total_num_check <- counts_dataframe |>
        mutate(row_number = row_number(),
               position_string = variantstring::position_from_variant_string(variant_string)) |>
        group_by(study_id, survey_id, position_string) |>
        summarise(total_num_check = sum(variant_num),
                  total_num = total_num[1],
                  row_number = row_number[1], .groups = "drop") |>
        mutate(overcount = total_num_check > total_num)
      
      if (any(total_num_check$overcount)) {
        
        # find problem rows
        problem_rows <- total_num_check |>
          filter(overcount == TRUE) |>
          pull(row_number)
        
        stop(paste0("in the counts_dataframe: for a given study, survey and variant, the sum of variant_num",
                    "\ncannot exceed the total_num. Problem rows in the counts_dataframe:\n",
                    paste(problem_rows, collapse = ", ")))
      }
      
    }
  )
)
