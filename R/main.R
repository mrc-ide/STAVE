
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
#' @import dplyr
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
    #' @param studies_dataframe a data.frame containing information at the study
    #'   level. This data.frame must have the following columns: study_id,
    #'   study_name, study_type, authors, publication_year, url
    #' @param surveys_dataframe a data.frame containing information at the
    #'   survey level. This data.frame must have the following columns:
    #'   study_key, survey_id, country_name, site_name, latitude, longitude, spatial_notes,
    #'   collection_start, collection_end, collection_day, time_notes.
    #' @param counts_dataframe a data.frame of genetic information. Must contain
    #'   the following columns: study_key, survey_key, variant_string, variant_num,
    #'   total_num.
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
      
      if (runif(1) < 0.05) {
        message("data correctly appended\nYou're smashing this :)")
      } else {
        message("data correctly appended")
      }
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
    #'   \code{TRUE}, which risks underestimating prevalence.
    #' @param return_full if \code{TRUE} (the default) returns the entire loaded
    #'   dataset, with prevalence equal to \code{NA} if there is no denominator.
    #'   If \code{FALSE} only returns entries for which there is a non-zero
    #'   denominator.
    #'   
    #'   @import dplyr
    get_prevalence = function(target_variant, keep_ambiguous = FALSE, prev_from_min = TRUE, return_full = TRUE) {
      
      # basic input checks
      assert_single_string(target_variant)
      assert_single_logical(keep_ambiguous)
      assert_single_logical(prev_from_min)
      
      # compare variant strings to get matches
      df_match <- variantstring::compare_variant_string(target_string = target_variant,
                                                        comparison_strings = private$counts$variant_string)
      
      df_prev <- private$counts |>
        mutate(match = df_match$match,
               ambiguous = df_match$ambiguous)
      
      # inform user if there are ambiguous matches to consider
      any_ambiguous <- any(df_prev$ambiguous)
      if ((keep_ambiguous == FALSE) & any_ambiguous) {
        warning(wrap_message(paste("There are some ambiguous matches to the target sequence due to mixed",
                                   "infections. Consider running with keep_ambiguous = TRUE to return both",
                                   "upper and lower bounds on the numerator. Currently ignoring all ambiguous",
                                   "matches, which may lead to prevalence being downwardly biased.",
                                   collapse = " ")))
      }
      
      # look for positional match, used in denominator calculation
      df_prev <- df_prev |>
        mutate(position_string = variantstring::position_from_variant_string(variant_string),
               match_pos = variantstring::compare_position_string(target_string = variantstring::position_from_variant_string(target_variant),
                                                                  comparison_strings = df_prev$variant_string))
      
      # calculate numerator and denominator of prevalence calculation
      df_prev_summary <- df_prev |>
        group_by(study_key, survey_key, position_string) |>
        summarise(numerator_min = sum(variant_num * match * (1 - ambiguous)),
                  numerator_max = sum(variant_num * match),
                  denominator = total_num[match_pos][1],
                  .groups = "drop") |>
        ungroup() |>
        mutate(denominator = ifelse(is.na(denominator), 0, denominator)) |>
        group_by(study_key, survey_key) |>
        summarise(numerator_min = sum(numerator_min),
                  numerator_max = sum(numerator_max),
                  denominator = sum(denominator),
                  .groups = "drop") |>
        ungroup()
      
      # calculate 95% exact binomial CI
      df_CI <- mapply(function(i) {
        if (df_prev_summary$denominator[i] == 0) {
          return(data.frame(prevalence = NA, prevalence_lower = NA, prevalence_upper = NA))
        } else {
          numerator <- ifelse(prev_from_min, df_prev_summary$numerator_min[i], df_prev_summary$numerator_max[i])
          denominator <- df_prev_summary$denominator[i]
          CI <- binom.test(x = numerator, n = denominator)$conf.int
          return(data.frame(prevalence = numerator / denominator* 1e2,
                            prevalence_lower = CI[1] * 1e2,
                            prevalence_upper = CI[2] * 1e2))
        }
      }, 1:nrow(df_prev_summary), SIMPLIFY = FALSE) |>
        bind_rows()
      
      df_prev_summary <- df_prev_summary |>
        bind_cols(df_CI)
      
      # tidy up names
      if (!keep_ambiguous) {
        df_prev_summary <- df_prev_summary |>
          rename(numerator = numerator_min) |>
          select(-numerator_max)
      }
      
      # optional trim to non-NA
      if (!return_full) {
        df_prev_summary <- df_prev_summary |>
          filter(denominator > 0)
      }
      
      # merge back with all study and survey info
      ret <- private$studies |>
        right_join(private$surveys |>
                     rename(study_id = study_key),
                   by = join_by(study_id)) |>
        right_join(df_prev_summary |>
                     rename(study_id = study_key,
                            survey_id = survey_key),
                   by = join_by(study_id, survey_id))
      
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
        select(study_key, survey_id) |>
        dplyr::filter(study_key %in% drop_study_id) |>
        rename(survey_key = survey_id)
      
      # drop IDs and save back into private objects
      private$studies <- self$get_studies() |>
        dplyr::filter(!(study_id %in% drop_study_id))
      
      private$surveys <- self$get_surveys() |>
        dplyr::filter(!(study_key %in% drop_study_id))
      
      private$counts <- self$get_counts() |>
        anti_join(drop_survey_id, by = join_by(study_key, survey_key))
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
      assert_ncol(studies_dataframe, 6)
      assert_eq(names(studies_dataframe), c("study_id", "study_name", "study_type", "authors", "publication_year", "url"))
      
      # study_id
      assert_non_NA(studies_dataframe$study_id)
      assert_valid_string(studies_dataframe$study_id, message_name = "study_id in studies_dataframe")
      
      # study_name
      w <- which(!is.na(studies_dataframe$study_name))
      if (any(w)) {
        assert_string(studies_dataframe$study_name[w])
      }
      
      # study_type
      assert_in(studies_dataframe$study_type, c("peer_reviewed", "preprint", "private", "other"))
      if (any(studies_dataframe$study_type == "private")) {
        warning("Some appended data are labelled as private. Be sure to remove these before making this data object public.")
      }
      
      # authors
      w <- which(!is.na(studies_dataframe$authors))
      if (any(w)) {
        assert_string(studies_dataframe$authors[w])
      }
      
      # publication_year
      w <- which(!is.na(studies_dataframe$publication_year))
      if (any(w)) {
        assert_pos_int(studies_dataframe$publication_year[w], zero_allowed = FALSE)
      }
      
      # url
      assert_non_NA(studies_dataframe$url)
      assert_string(studies_dataframe$url)
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
      assert_ncol(surveys_dataframe, 11)
      assert_eq(names(surveys_dataframe), c("study_key", "survey_id", "country_name", "site_name", "latitude",
                                            "longitude", "spatial_notes", "collection_start", "collection_end",
                                            "collection_day", "time_notes"))
      
      # study_key
      assert_non_NA(surveys_dataframe$study_key)
      assert_valid_string(surveys_dataframe$study_key, message_name = "study_key in surveys_dataframe")
      
      # survey_id
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
      
      # latitude
      assert_non_NA(surveys_dataframe$latitude)
      assert_bounded(surveys_dataframe$latitude, left = -180, right = 180)
      
      # longitude
      assert_non_NA(surveys_dataframe$longitude)
      assert_bounded(surveys_dataframe$longitude, left = -180, right = 180)
      
      # spatial_notes
      w <- which(!is.na(surveys_dataframe$spatial_notes))
      if (any(w)) {
        assert_string(surveys_dataframe$spatial_notes[w])
      }
      
      # collection_start
      w <- which(!is.na(surveys_dataframe$collection_start))
      if (any(w)) {
        assert_valid_ymd(surveys_dataframe$collection_start[w])
      }
      
      # collection_end
      w <- which(!is.na(surveys_dataframe$collection_end))
      if (any(w)) {
        assert_valid_ymd(surveys_dataframe$collection_end[w])
      }
      
      # collection_day
      assert_non_NA(surveys_dataframe$collection_day)
      assert_valid_ymd(surveys_dataframe$collection_day)
      
      # time_notes
      w <- which(!is.na(surveys_dataframe$time_notes))
      if (any(w)) {
        assert_string(surveys_dataframe$time_notes[w])
      }
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
      assert_eq(names(counts_dataframe), c("study_key", "survey_key", "variant_string", "variant_num", "total_num"))
      
      # study_key
      assert_non_NA(counts_dataframe$study_key)
      assert_valid_string(counts_dataframe$study_key, message_name = "study_key in counts_dataframe")
      
      # survey_key
      assert_non_NA(counts_dataframe$survey_key)
      assert_valid_string(counts_dataframe$survey_key, message_name = "survey_key in counts_dataframe")
      
      # variant_string
      assert_non_NA(counts_dataframe$variant_string)
      variantstring::check_variant_string(counts_dataframe$variant_string)
      
      # variant_num
      assert_non_NA(counts_dataframe$variant_num)
      assert_pos_int(counts_dataframe$variant_num, zero_allowed = TRUE)
      
      # total_num
      assert_non_NA(counts_dataframe$total_num)
      assert_pos_int(counts_dataframe$total_num, zero_allowed = FALSE)
      
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
      # - survey_ids are unique within each study (but not necessarily completely unique)
      # - study_keys can be found in studies table
      # - combination of study_ids and survey_ids are referenced at least once in counts table
      
      # counts_dataframe:
      # - combination of study_key and survey_key can be found in surveys table
      # - variant_strings are unique within each study-survey (even if listed in different order)
      # - total_num is consistent between variants that correspond to the same positions
      # - after splitting by study-survey and also by gene-locus combination (i.e. the first two characteristics
      #   of each encoded gene, but not by the observed amino acids):
      #       - variant_num cannot sum to more than total_num
      
      # studies_dataframe:
      assert_noduplicates(studies_dataframe$study_id)
      if (any(studies_dataframe$study_id %in% private$studies$study_id)) {
        stop(wrap_message(paste("studies_dataframe cannot contain any study_ids that already exist in the loaded",
                                "studies table. If you want to add to an existing study then you must first drop the",
                                "currently loaded version. See ?drop_study() for how to do this", collapse = "")))
      }
      if (!all(studies_dataframe$study_id %in% surveys_dataframe$study_key)) {
        stop(wrap_message(paste("every study_id in the studies_dataframe must be referenced at least once in the",
                                "study_key column of the surveys_dataframe", collapse = "")))
      }
      
      # surveys_dataframe:
      # make unique ID as combination of study and survey IDs
      surveys_dataframe <- surveys_dataframe |>
        mutate(combined_ID = paste(study_key, survey_id, sep = ":"))
      
      assert_noduplicates(surveys_dataframe$combined_ID, message = wrap_message(paste("survey_ids are allowed to be re-used between studies,",
                                                                                      "but the same survey_id cannot be used within a study",
                                                                                      collapse = "")))
      if (!all(surveys_dataframe$study_key %in% studies_dataframe$study_id)) {
        stop(wrap_message(paste("every study_key in the surveys_dataframe must be present as a study_id in the",
                                "studies_dataframe", collapse = "")))
      }
      if (!all(surveys_dataframe$survey_id %in% counts_dataframe$survey_key)) {
        stop(wrap_message(paste("every survey_id in the surveys_dataframe must be referenced at least once in the",
                                "survey_key column of the counts_dataframe", collapse = "")))
      }
      
      # counts_dataframe:
      # make unique ID as combination of study and survey IDs
      counts_dataframe <- counts_dataframe |>
        mutate(combined_ID = paste(study_key, survey_key, sep = ":"))
      
      if (!all(counts_dataframe$combined_ID %in% surveys_dataframe$combined_ID)) {
        stop(wrap_message(paste("every study_key-survey_key combination in the counts_dataframe must be present",
                                "as a combination in the surveys_dataframe", collapse = "")))
      }
      
      # ensure no variant is present more than once (should not be after sorting alphabetically)
      variant_duplicated <- counts_dataframe |>
        mutate(variant_string_order = variantstring::order_variant_string(variant_string)) |>
        group_by(combined_ID, variant_string_order) |>
        summarise(variant_duplicated = any(duplicated(variant_string_order)), .groups = "drop") |>
        pull(variant_duplicated)
      
      if (any(variant_duplicated)) {
        stop(wrap_message(paste("the exact same variant cannot be present more than once in the same survey.",
                                "This includes the same variant with the genes listed in different order", collapse = "")))
      }
      
      # ensure total_num is consistent between variants that correspond to the same positions
      distinct_total_num <- counts_dataframe |>
        dplyr::mutate(position_string = variantstring::position_from_variant_string(variant_string)) |>
        dplyr::group_by(combined_ID, position_string) |>
        dplyr::summarise(distinct_total_num = length(unique(total_num)), .groups = "drop") |>
        dplyr::pull(distinct_total_num)
      
      if (any(distinct_total_num > 1)) {
        stop(wrap_message(paste("if there are multiple rows in a survey that correspond to the same gene-locus",
                                "combination then they must have the same total_num. This includes variants",
                                "with heterozygous calls", collapse = "")))
      }
      
      # sum variants and determine if any where counts exceed the total_num
      total_num_check <- counts_dataframe |>
        dplyr::mutate(row_number = dplyr::row_number(),
                      position_string = variantstring::position_from_variant_string(variant_string)) |>
        dplyr::group_by(combined_ID, position_string) |>
        dplyr::summarise(total_num_check = sum(variant_num),
                         total_num = total_num[1],
                         row_number = row_number[1], .groups = "drop") |>
        dplyr::mutate(overcount = total_num_check > total_num)
      
      if (any(total_num_check$overcount)) {
        
        # find problem rows
        problem_rows <- total_num_check |>
          dplyr::filter(overcount == TRUE) |>
          dplyr::pull(row_number)
        
        message(paste("in the counts_dataframe: for a given study, survey and variant, the sum of variant_num",
                      "\ncannot exceed the total_num. Problem rows in the counts_dataframe:\n",
                      paste(problem_rows, collapse = ", "), collapse = ""))
        
        stop()
      }
      
    }
  )
)


