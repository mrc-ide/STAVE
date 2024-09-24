
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
#'     and author names. Each study is indexed with a unique study_ID.
#'     \item surveys: Information on the surveys represented within a study. A
#'     survey is defined here as a discrete instance of data collection, which
#'     includes information on geography (latitude and longitude) and collection
#'     time. Surveys are given survey_IDs and are linked to a particular study
#'     through the study_ID.
#'     \item counts: The actual genetic information, which is linked to a
#'     particular survey through the survey_ID. Genetic variants are encoded in
#'     character strings that must follow a specified format, and the number of
#'     times this variant was observed among the total sample is stored in
#'     columns.
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
        cat(sprintf("Studies: %s", nrow(private$studies)))
        n_surveys <- length(unique(private$surveys$survey_ID))
        cat(sprintf("\nSurveys: %s", n_surveys))
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
    #' Append new data
    #' @param studies_dataframe a data.frame containing information at the study
    #'   level. This data.frame must have the following columns: study_ID,
    #'   study_name, study_type, authors, publication_year, url
    #' @param surveys_dataframe a data.frame containing information at the
    #'   survey level. This data.frame must have the following columns:
    #'   study_key, survey_ID, country_name, site_name, lat, lon, spatial_notes,
    #'   collection_start, collection_end, time_notes. The study_key element
    #'   must correspond to a study_ID in the studies_dataframe.
    #' @param counts_dataframe a data.frame of genetic information. Must contain
    #'   the following columns: survey_key, variant_string, variant_num,
    #'   total_num. The survey_key element must correspond to a valid survey_ID
    #'   in the surveys_dataframe.
    append_data = function(studies_dataframe, surveys_dataframe, counts_dataframe) {
      
      # basic checks on format of input data.frames. These checks ensure that
      # each independent row of the input data are correctly formatted
      private$check_format_studies(studies_dataframe)
      private$check_format_surveys(surveys_dataframe)
      private$check_format_counts(counts_dataframe)
      private$check_variant_string(counts_dataframe$variant_string)
      
      # structural checks on data. These are more focused on the relationships
      # between data, for example ensuring no duplicate study_IDs, or ensuring
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
    #' @param target_variant the name of the variant on which we want to
    #'   calculate prevalence, for example crt:72:C. Note that there can be no
    #'   heterozygous calls within this name.
    #' @param keep_ambiguous there may be variants in the data for which the
    #'   target_variant could be in the sample, but this cannot be proven
    #'   conclusively. For example, the sequence A_A_A may be a match to the
    #'   sequence A/C_A/C_A or it may not, these are unphased genotypes so we
    #'   cannot be sure. If \code{keep_ambiguous = TRUE} then both a min and a
    #'   max numerator are reported that either exclude all ambiguous calls
    #'   (min) or include all ambiguous calls (max). If \code{FALSE} (the
    #'   default) then only the min is reported.
    #' @param prev_from_min the output object includes a point estimate of the
    #'   prevalence along with exact binomial confidence intervals. These must
    #'   be calculated from one of \code{numerator_min} or \code{numerator_max}
    #'   in the case of ambiguous calls. This argument sets which one of these
    #'   numerators is used in the calculation.
    #'   
    #'   @import dplyr
    get_prevalence = function(target_variant, keep_ambiguous = FALSE, prev_from_min = TRUE) {
      
      # basic input checks
      assert_single_string(target_variant)
      assert_single_logical(keep_ambiguous)
      assert_single_logical(prev_from_min)
      
      # check that the target_variant is in long-string format
      private$check_variant_string(target_variant)
      
      # can only calculate prev from max if keeping ambiguous
      if (!keep_ambiguous & !prev_from_min) {
        stop(wrap_message("keep_ambiguous must be TRUE in order to calculate from max values (i.e. when using prev_from_min = FALSE)"))
      }
      
      # check no hets in target_variant
      if (grepl("[/|]", target_variant)) {
        stop("No het calls are allowed in the target sequence")
      }
      
      # get version of target_variant without the amino acid sequence
      target_variant_drop <- drop_amino(target_variant)
      
      # combine all three data.frames into one using IDs and keys
      counts_ID <- dplyr::rename(private$counts, survey_ID = survey_key)
      surveys_ID <- dplyr::rename(private$surveys, study_ID = study_key)
      df_combined <- private$studies |>
        dplyr::right_join(surveys_ID, by = "study_ID") |>
        dplyr::right_join(counts_ID, by = "survey_ID")
      
      # get stripped variant from the data (minus amino acid) as we will need to
      # group by this variable. For example, crt_1:A would match with both
      # crt:1:A and crt:1:A/C, and numerator (but not denominator) needs to be
      # summed over both of these.
      df_combined <- df_combined |>
        dplyr::mutate(variant_drop = drop_amino(counts_ID$variant_string, sort_gene_name = TRUE))
      
      # find if the target variant matches including the amino acid sequence
      # (for the numerator) and without the amino acid sequence (for the
      # denominator). Note that numerator matches may include "Yes", "No" or
      # "Ambiguous"
      df_combined <- df_combined |>
        dplyr::mutate(variant_match = private$get_variant_matches(target_variant = target_variant),
                      variant_drop_match = private$get_variant_matches(target_variant = target_variant_drop))
      
      # inform user if there are ambiguous counts to consider
      any_ambiguous <- any(df_combined$variant_match == "Ambiguous")
      if ((keep_ambiguous == FALSE) & any_ambiguous) {
        warning(wrap_message(paste("There are some ambiguous matches to the target sequence due to mixed",
                                   "infections. Consider running with keep_ambiguous = TRUE to return both",
                                   "upper and lower bounds on the numerator. Currently ignoring all ambiguous",
                                   "matches, which may lead to prevalence being biased downward",
                                   collapse = " ")))
      }
      
      # perform filtering, grouping, and first step aggregation of numerator
      # (not yet summing variants at the same locus)
      # This version is for unambiguous calls only.
      df_numerator_unambiguous_stage1 <- df_combined |>
        dplyr::mutate(row_number = dplyr::row_number()) |>
        dplyr::filter(variant_match == "Yes") |>
        dplyr::group_by(study_ID, survey_ID, variant_drop) |>
        dplyr::summarise(numerator_min = sum(variant_num),
                         row_number = row_number[1],
                         .groups = "drop") |>
        dplyr::ungroup()
      
      # second stage numerator aggregation where we sum over variants at the
      # same locus
      df_numerator_unambiguous <- df_numerator_unambiguous_stage1 |>
        dplyr::group_by(study_ID, survey_ID) |>
        dplyr::summarise(numerator_min = sum(numerator_min),
                         row_number = row_number[1],
                         .groups = "drop") |>
        dplyr::ungroup() |>
        dplyr::arrange(row_number) |>
        dplyr::select(-row_number)
      
      # perform filtering, grouping, and first step aggregation of denominator
      # (not yet summing variants at the same locus)
      df_denominator_stage1 <- df_combined |>
        dplyr::mutate(row_number = dplyr::row_number()) |>
        dplyr::filter(variant_drop_match == "Yes") |>
        dplyr::group_by(study_ID, survey_ID, variant_drop) |>
        dplyr::summarise(denominator = total_num[1],
                         row_number = row_number[1],
                         .groups = "drop") |>
        dplyr::ungroup()
      
      # second stage denominator aggregation where we sum over variants at the
      # same locus
      df_denominator <- df_denominator_stage1 |>
        dplyr::group_by(study_ID, survey_ID) |>
        dplyr::summarise(denominator = sum(denominator),
                         row_number = row_number[1],
                         .groups = "drop") |>
        dplyr::ungroup() |>
        dplyr::arrange(row_number) |>
        dplyr::select(-row_number)
      
      # copy over object to df_numerator. The ambiguous counts will be
      # optionally merged with this object
      df_numerator <- df_numerator_unambiguous
      
      # optionally do the same for ambiguous calls
      if (keep_ambiguous) {
        
        # repeat the process of aggregation for ambiguous calls
        df_numerator_ambiguous_stage1 <- df_combined |>
          dplyr::mutate(row_number = dplyr::row_number()) |>
          dplyr::filter(variant_match  %in% c("Yes", "Ambiguous")) |>
          dplyr::group_by(study_ID, survey_ID, variant_drop) |>
          dplyr::summarise(numerator_max = sum(variant_num),
                           row_number = row_number[1],
                           .groups = "drop") |>
          dplyr::ungroup()
        
        # second stage numerator aggregation where we sum over variants at the
        # same locus
        df_numerator_ambiguous <- df_numerator_ambiguous_stage1 |>
          dplyr::group_by(study_ID, survey_ID) |>
          dplyr::summarise(numerator_max = sum(numerator_max),
                           row_number = row_number[1],
                           .groups = "drop") |>
          dplyr::ungroup() |>
          dplyr::arrange(row_number) |>
          dplyr::select(-row_number)
        
        # merge with unambiguous
        df_numerator <- df_numerator |>
          dplyr::left_join(df_numerator_ambiguous, by = dplyr::join_by(study_ID, survey_ID))
      }
      
      # merge numerator with denominator
      df_prev <- df_numerator |>
        dplyr::left_join(df_denominator, by = dplyr::join_by(study_ID, survey_ID))
      
      # calculate prevalence point estimate and 95% exact binomial confidence
      # intervals
      if (nrow(df_prev) == 0) {
        df_prev <- df_prev |>
          dplyr::mutate(prevalence = NA,
                        prevalence_lower = NA,
                        prevalence_upper = NA)
      } else {
        if (prev_from_min) {
          prev_CI <- mapply(function(x, n) binom.test(x, n)$conf.int, df_prev$numerator_min, df_prev$denominator)
          df_prev <- df_prev |>
            dplyr::mutate(prevalence = numerator_min / denominator,
                          prevalence_lower = prev_CI[1,],
                          prevalence_upper = prev_CI[2,])
          
        } else {
          prev_CI <- mapply(function(x, n) binom.test(x, n)$conf.int, df_prev$numerator_max, df_prev$denominator)
          df_prev <- df_prev |>
            dplyr::mutate(prevalence = numerator_max / denominator,
                          prevalence_lower = prev_CI[1,],
                          prevalence_upper = prev_CI[2,])
        }
      }
      
      # tidy up names
      if (!keep_ambiguous) {
        df_prev <- df_prev |>
          dplyr::rename(numerator = numerator_min)
      }
      
      # merge back with all study and survey info
      ret <- private$studies |>
        dplyr::right_join(surveys_ID, by = "study_ID") |>
        dplyr::left_join(df_prev, by = dplyr::join_by(study_ID, survey_ID))
      
      # replace some NAs with zeros
      if (keep_ambiguous) {
        ret$numerator_min <- ifelse(is.na(ret$numerator_min), 0, ret$numerator_min)
        ret$numerator_max <- ifelse(is.na(ret$numerator_max), 0, ret$numerator_max)
      } else {
        ret$numerator <- ifelse(is.na(ret$numerator), 0, ret$numerator)
      }
      ret$denominator <- ifelse(is.na(ret$denominator), 0, ret$denominator)
      
      return(ret)
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
    # check basic formatting (variables types etc) of studies data.frame
    check_format_studies = function(studies_dataframe) {
      
      # CHECKS
      # - is data.frame
      # - has exactly the right number of columns and correctly named
      # - variables have the correct type and range
      # - check that character string variables are valid
      
      assert_dataframe(studies_dataframe)
      assert_ncol(studies_dataframe, 6)
      assert_eq(names(studies_dataframe), c("study_ID", "study_name", "study_type", "authors", "publication_year", "url"))
      assert_string(studies_dataframe$study_ID)
      assert_valid_string(studies_dataframe$study_ID, message_name = "study_ID in studies_dataframe")
      if (!all(is.na(studies_dataframe$study_name))) {
        assert_string(studies_dataframe$study_name)
      }
      assert_in(studies_dataframe$study_type, c("peer_reviewed", "preprint", "other"))
      if (!all(is.na(studies_dataframe$authors))) {
        assert_string(studies_dataframe$authors)
      }
      if (!all(is.na(studies_dataframe$publication_year))) {
        assert_pos_int(studies_dataframe$publication_year, zero_allowed = FALSE)
      }
      assert_string(studies_dataframe$url)
    },
    
    # -----------------------------------
    # check basic formatting (variables types etc) of surveys data.frame
    check_format_surveys = function(surveys_dataframe) {
      
      # CHECKS
      # - is data.frame
      # - has exactly the right number of columns and correctly named
      # - variables have the correct type and range. Allow for missing data in some cases
      # - check that character string variables are valid
      # - country_name is within the valid list of names
      # - collection_start must be before collection_end
      # - collection_day is between collection_start and collection_end
      
      assert_dataframe(surveys_dataframe)
      assert_ncol(surveys_dataframe, 11)
      assert_eq(names(surveys_dataframe), c("study_key", "survey_ID", "country_name", "site_name", "lat", "lon", "spatial_notes",
                                            "collection_start", "collection_end", "collection_day", "time_notes"))
      assert_string(surveys_dataframe$survey_ID)
      assert_valid_string(surveys_dataframe$survey_ID, message_name = "survey_ID in surveys_dataframe")
      assert_string(surveys_dataframe$study_key)
      assert_valid_string(surveys_dataframe$study_key, message_name = "study_key in surveys_dataframe")
      
      if (!all(is.na(surveys_dataframe$country_name))) {
        assert_string(surveys_dataframe$country_name)
      }
      if (!all(is.na(surveys_dataframe$site_name))) {
        assert_string(surveys_dataframe$site_name)
      }
      assert_bounded(surveys_dataframe$lat, left = -180, right = 180)
      assert_bounded(surveys_dataframe$lon, left = -180, right = 180)
      if (!all(is.na(surveys_dataframe$spatial_notes))) {
        assert_string(surveys_dataframe$spatial_notes)
      }
      if (!all(is.na(surveys_dataframe$collection_start))) {
        assert_string(surveys_dataframe$collection_start)
      }
      if (!all(is.na(surveys_dataframe$collection_end))) {
        assert_string(surveys_dataframe$collection_end)
      }
      assert_valid_ymd(surveys_dataframe$collection_day)
      if (!all(is.na(surveys_dataframe$time_notes))) {
        assert_string(surveys_dataframe$time_notes)
      }
    },
    
    # -----------------------------------
    # check basic formatting (variables types etc) of counts data.frame
    check_format_counts = function(counts_dataframe) {
      
      # CHECKS
      # - is data.frame
      # - has exactly the right number of columns and correctly named
      # - variables have the correct type and range. Allow for missing data in some cases
      # - check that character string variables are valid. Note that a more thorough series of
      #   checks on the variant_string format will occur separately through the
      #   check_variant_string() function
      # - variant_num cannot exceed total_num
      
      assert_dataframe(counts_dataframe)
      assert_ncol(counts_dataframe, 4)
      assert_eq(names(counts_dataframe), c("survey_key", "variant_string", "variant_num", "total_num"))
      assert_string(counts_dataframe$survey_key)
      assert_valid_string(counts_dataframe$survey_key, message_name = "survey_key in counts_dataframe")
      assert_string(counts_dataframe$variant_string)
      assert_pos_int(counts_dataframe$variant_num, zero_allowed = TRUE)
      assert_pos_int(counts_dataframe$total_num, zero_allowed = FALSE)
      assert_leq(counts_dataframe$variant_num, counts_dataframe$total_num, message = "variant_num cannot exceed total_num")
      
    },
    
    # -----------------------------------
    # check that the variant_string is correctly formatted. This is a complex
    # character string that must obey strict rules
    check_variant_string = function(variant_string) {
      
      # CHECKS
      # - after splitting by gene (semicolon), and then characteristic (colon), contains exactly three characteristics per gene
      # - gene name (first characteristic):
      #     - contains only allowed characters
      # - codon position (second characteristic):
      #     - contains only allowed characters
      #     - split by codon (underscore) and count number of codons. Will use this number in later check
      #     - codons are positive integers and sorted increasing
      # - amino acid (third characteristic):
      #     - contains only allowed characters
      #     - remove underscores
      #     - split into groups, allowing for / and | symbols
      #     - number of groups matches number of codons
      #     - no duplicated adjacent / or | symbols. No symbols at start or end of group
      #     - a group cannot contain a mix of both phased and unphased amino acids (| and /)
      #     - all groups with phased heterozygotes contain the same number of amino acids. This is true
      #       across phased genes. Mixed phased and unphased hets are allowed (| and /), in which case
      #       unphased hets are not required to have the same number of amino acids as phased hets.
      #     - hets cannot have the same amino acid multiple times
      # - no duplicate gene names
      
      # create objects for storing which rows fail. This allows for more
      # informative error messages in which several issues can be listed at once
      # rather than breaking at the first issue
      n_string <- length(variant_string)
      valid <- rep(TRUE, n_string)
      reason <- rep(NA, n_string)
      
      # get list of valid amino acid characters. Do this once here to avoid
      # repetition for each row
      IUPAC_df <- allowed_amino_acids()
      valid_amino_characters <- paste0("^[", paste(IUPAC_df$IUPAC_amino_acid_code, collapse = ""), "_/|]+$")
      
      # input may be a single string or a vector of strings. Loop through elements
      for (i in 1:n_string) {
        
        # split into genes
        s1 <- stringr::str_split(variant_string[i], ";", n = Inf, simplify = FALSE)[[1]]
        n_genes <- length(s1)
        #s1
        
        if (any(s1 == "")) {
          valid[i] <- FALSE
          reason[i] <-"contains one more empty genes"
          next()
        }
        
        # keep track of gene names
        gene_names <- rep(NA, n_genes)
        
        # keep track of the number of distinct amino acid calls in phased
        # heterozygous loci. This must be the same over all genes
        gene_n_phased <- rep(NA, n_genes)
        
        # loop over genes
        for (j in seq_along(s1)) {
          
          # split this gene into characteristics
          x <- stringr::str_split(s1[j], ":", n = Inf, simplify = FALSE)[[1]]
          gene_names[j] <- x[1]
          #x
          
          if (length(x) != 3) {
            valid[i] <- FALSE
            reason[i] <- sprintf("gene %s does not contain three characteristics separated by a colon", j)
            next()
          }
          
          if (!grepl("^[a-z][a-z0-9]*$", x[1])) {
            valid[i] <- FALSE
            reason[i] <- sprintf("gene name %s contains invalid characters", x[1])
            next()
          }
          
          if (!grepl("^[0-9_]+$", x[2])) {
            valid[i] <- FALSE
            reason[i] <- "codon positions contain invalid characters"
            next()
          }
          
          # split second characteristic by codon and make numeric
          y <- stringr::str_split(x[2], "_", n = Inf, simplify = FALSE)[[1]] |>
            as.numeric()
          #y
          
          if (any(y == 0)) {
            valid[i] <- FALSE
            reason[i] <- "codon position cannot be zero"
            next()
          }
          
          if (any(duplicated(y))) {
            valid[i] <- FALSE
            reason[i] <- "codon contains duplicated positions"
            next()
          }
          
          if (!all(diff(y) > 0)) {
            valid[i] <- FALSE
            reason[i] <- "codon positions must be sorted increasing"
            next()
          }
          
          if (!grepl(valid_amino_characters, x[3])) {
            valid[i] <- FALSE
            reason[i] <- "amino acid sequence contains invalid characters. See ?allowed_amino_acids()"
            next()
          }
          
          # remove underscores
          no_underscore <- gsub("_", "", x[3])
          #no_underscore
          
          if (grepl("^[/|]|[/|]$", no_underscore)) {
            valid[i] <- FALSE
            reason[i] <- "amino acid sequence cannot start or end with / or | symbols"
            next()
          }
          
          if (grepl("[/|]{2,}", no_underscore)) {
            valid[i] <- FALSE
            reason[i] <- "amino acid sequence cannot contain adjacent / or | symbols"
            next()
          }
          
          # split third characteristic using a regular expression to group
          # sequences of letters connected by '/' or '|'
          z <- unlist(strsplit(no_underscore, "(?<=[A-Z])(?=[^/|])", perl = TRUE))
          #z
          
          if (length(z) != length(y)) {
            valid[i] <- FALSE
            reason[i] <- sprintf("number of amino acid loci (%s) must equal the number of codon positions (%s)",
                                 length(z), length(y))
            next()
          }
          
          if (any(grepl("/", z) & grepl("\\|", z))) {
            valid[i] <- FALSE
            reason[i] <- paste("a single amino acid locus cannot contain both the / and",
                               "the | symbols. Heterozygous loci must be either entirely phased or",
                               "entirely unphased", collapse = "")
            next()
          }
          
          # checks on phased heterozygous loci (if present)
          z_phased_het <- z[grepl("\\|", z)]
          if (length(z_phased_het) != 0) {
            
            n_phased <- mapply(length, strsplit(z_phased_het, "\\|", ))
            if (length(unique(n_phased)) != 1) {
              valid[i] <- FALSE
              reason[i] <- paste("if there are phased heterozygous loci within a gene then the number of",
                                 "distinct amino acids must be the same over all of these loci. Variable",
                                 "degrees of phasing are not allowed", collapse = "")
              next()
            }
            
            # store number of distinct amino acids in phased calls for this
            # gene. This will be tested for consistency across genes
            gene_n_phased[j] <- n_phased[1]
          }
          
          # split into distinct amino acids
          q <- strsplit(z, "[/|]", )
          #q
          
          if (any(mapply(function(x) any(duplicated(x)), q))) {
            valid[i] <- FALSE
            reason[i] <- "the same amino acid cannot be present more than once in a heterozygous call"
            next()
          }
          
        } # end loop over genes
        
        if (length(unique(na.omit(gene_n_phased))) > 1) {
          valid[i] <- FALSE
          reason[i] <- paste("if there are phased heterozygous loci then the same number of distinct amino",
                             "acids must be present all of these loci AND over all genes", collapse = "")
          next()
        }
        
        if (any(duplicated(gene_names))) {
          valid[i] <- FALSE
          reason[i] <- "duplicated gene names"
          next()
        }
        
      } # end loop over elements of variant_string
      
      # if any invalid entries then print all issues before exit
      if (any(!valid)) {
        message("The following issues were found in variant_string:")
        for (i in seq_along(valid)) {
          if (!valid[i]) {
            message(wrap_message(sprintf("  - row %s: %s", i, reason[i])))
          }
        }
        
        stop()
      }
      
    },  # end check_variant_string()
    
    # -----------------------------------
    # check that new inputs can be appended without breaking overall structure
    check_structure = function(studies_dataframe, surveys_dataframe, counts_dataframe) {
      
      # CHECKS
      # studies_dataframe:
      # - new study_IDs are unique
      # - new study_IDs are not already present in old study_IDs
      # - new study_IDs are referenced at least once in new study_keys
      
      # surveys_dataframe:
      # - new survey_IDs are unique within each study
      # - new study_keys map to real new study_IDs
      # - new survey_IDs are referenced at least once in survey_keys within new counts table
      
      # counts_dataframe:
      # - new survey_keys map to real new survey_IDs
      # - no duplicate variants after split by study_ID
      # - split by survey_key and also by gene-locus combination (i.e. the first two characteristics
      #   of each encoded gene, but not by the observed amino acids):
      #       - total_num must be identical over rows
      #       - variant_num cannot sum to more than total_num
      
      # studies_dataframe:
      assert_noduplicates(studies_dataframe$study_ID)
      if (any(studies_dataframe$study_ID %in% private$studies$study_ID)) {
        stop(wrap_message(paste("studies_dataframe cannot contain any study_IDs that already exist in the loaded",
                                "studies table. If you want to add to an existing study then you must first drop the",
                                "currently loaded version. See ?drop_study() for how to do this", collapse = "")))
      }
      if (!all(studies_dataframe$study_ID %in% surveys_dataframe$study_key)) {
        stop(wrap_message(paste("every study_ID in the studies_dataframe must be referenced at least once in the",
                                "study_key column of the surveys_dataframe", collapse = "")))
      }
      
      # surveys_dataframe:
      # make unique ID as combination of study and survey IDs
      UID <- paste(surveys_dataframe$study_key, surveys_dataframe$survey_ID)
      assert_noduplicates(UID, message = wrap_message(paste("survey_IDs are allowed to be re-used between studies,",
                                                            "but the same survey_ID cannot be used within a study",
                                                            collapse = "")))
      if (!all(surveys_dataframe$study_key %in% studies_dataframe$study_ID)) {
        stop(wrap_message(paste("every study_key in the surveys_dataframe must be present as a study_ID in the",
                                "studies_dataframe", collapse = "")))
      }
      if (!all(surveys_dataframe$survey_ID %in% counts_dataframe$survey_key)) {
        stop(wrap_message(paste("every survey_ID in the surveys_dataframe must be referenced at least once in the",
                                "survey_key column of the counts_dataframe", collapse = "")))
      }
      
      # counts_dataframe:
      if (!all(counts_dataframe$survey_key %in% surveys_dataframe$survey_ID)) {
        stop(wrap_message(paste("every survey_key in the counts_dataframe must be present as a survey_ID in the",
                                "surveys_dataframe", collapse = "")))
      }
      
      if (any(duplicated(sort_gene_name(counts_dataframe$variant_string)))) {
        stop(wrap_message(paste("the exact same variant cannot be present more than once in the same survey.",
                                "This includes the same variant with the genes listed in different order", collapse = "")))
      }
      
      # get number of distinct variants per survey
      distinct_variants <- counts_dataframe |>
        dplyr::mutate(variant_drop = drop_amino(counts_dataframe$variant_string)) |>
        dplyr::group_by(survey_key, variant_drop) |>
        dplyr::summarise(distinct_total_num = length(unique(total_num)), .groups = "drop") |>
        dplyr::pull(distinct_total_num)
      
      if (any(distinct_variants > 1)) {
        stop(wrap_message(paste("If there are multiple rows in a survey that correspond to the same gene-locus",
                                "combination then they must have the same total_num. This includes variants",
                                "with heterozygous calls", collapse = "")))
      }
      
      # sum variants and determine if any where counts exceed the total_num
      total_num_check <- counts_dataframe |>
        dplyr::mutate(variant_drop = drop_amino(counts_dataframe$variant_string)) |>
        dplyr::group_by(survey_key, variant_drop) |>
        dplyr::summarise(total_num_check = sum(variant_num),
                         total_num = total_num[1], .groups = "drop") |>
        dplyr::mutate(overcount = total_num_check > total_num) |>
        dplyr::pull(overcount)
      
      if (any(total_num_check)) {
        stop("For a given survey and variant, the sum of variant_num cannot exceed the total_num")
      }
      
    },
    
    # -----------------------------------
    # return a vector of whether the target_variant matches each entry in
    # counts$variant_string. Matches can be "Yes", "No", or "Ambiguous"
    get_variant_matches = function(target_variant) {
      
      # get target variant into list format
      vt <- variant_to_list(target_variant)
      
      # store whether each row in data matches the target variant. Matches can
      # be "Yes", "No", or "Ambiguous"
      match_vec <- rep("No", length(private$counts$variant_string))
      
      # loop through genes
      for (i in seq_along(private$counts$variant_string)) {
        
        # get this data variant into list format
        vd <- variant_to_list(private$counts$variant_string[i])
        
        # all target gene names must be in data
        if (!all(names(vt) %in% names(vd))) {
          next()
        }
        
        # loop through target genes
        match_gene <- rep("No", length(vt))
        for (j in seq_along(vt)) {
          
          # find matching gene
          w <- which(names(vd) == names(vt)[j])
          
          # all target positions must be in data
          if (!all(vt[[j]]$pos %in% vd[[w]]$pos)) {
            next()
          }
          
          # if not matching on amino acid sequence then criteria already met
          if (is.null(vt[[j]]$amino)) {
            match_gene[j] <- "Yes"
            next()
          }
          
          # extract target and data amino acids for convenience
          m <- match(vt[[j]]$pos, vd[[w]]$pos)
          at <- vt[[j]]$amino
          ad <- vd[[w]]$amino[m]
          
          # matching on phased and unphased heterozygous loci is where this gets
          # tricky. Matches can be 1) "No", 2) "Yes" if the target sequence is
          # definitely within the data, or 3) "Ambiguous" if the target sequence
          # may be within the data but we cannot know for sure.
          # 
          # if there is a single unphased het call within the data then
          # we know that both haplotypes were present. Therefore we should match
          # on either of them. For example, A_A_A should match A_A/C_A because
          # we know the target haplotype must be present within this mixed
          # infection. Therefore, this is either "Yes" or "No" for a match.
          # 
          # If there is more than one unphased het call this is no longer the
          # case. Assume COI=2 and there are 2 unphased het loci. Then there are
          # 2 possible combinations of haplotypes that cannot be distinguished.
          # For example, A_A/C_A/C could be caused by A_A_A + A_C_C, or equally
          # by A_A_C + A_C_A. If the target is A_A_A then this is only a match
          # in the first case and not the second. If COI=4 then this problem
          # would go away because, once again, we would be certain of a match.
          # But in the absence of COI information we must be pessimistic and
          # assume the worst case, which is that COI=2 and we cannot distinguish
          # haplotypes. Hence, A_A_A should be counted as an "Ambiguous" match
          # against A_A/C_A/C.
          # 
          # Moving on to phased het calls, the number of haplotypes in the
          # mixture is simply equal to the degree of the het call (i.e. the
          # number of amino acids). This is true for any number of het sites;
          # unlike for unphased het calls it does not get more uncertain with
          # more loci. For example, A_A|C_A|C definitely contains the two
          # haplotypes A_A_A and A_C_C. It should therefore be a definite match
          # against either of these target sequences.
          #
          # However, if we were to add in a single unphased locus the
          # combinations would expand. For example, A\C_A|C_A|C could be caused
          # by A_A_A + C_C_C or C_A_A + A_C_C. Therefore A_A_A should be counted
          # as an ambiguous match against this sequence.
          # 
          # Bringing this together:
          #   - if the sequence contains no het calls then this allows for an
          #     unambiguous match.
          #   - if the sequence contains one unphased het call but no more het
          #     calls (whether phased or unphased) then this allows for an
          #     unambiguous match.
          #   - if the sequence contains only phased het calls then this allows
          #     for an unambiguous match.
          #   - in all other situations, this allows for an ambiguous match
          
          # is this a match
          is_match <- all(mapply(function(at, ad) grepl(at, ad), at, ad))
          
          # is this ambiguous
          num_unphased_het <- sum(grepl("/", ad))
          num_phased_het <- sum(grepl("\\|", ad))
          is_ambiguous <- (num_unphased_het > 1) || ((num_unphased_het > 0) && (num_phased_het > 0))
          
          if (is_match) {
            if (is_ambiguous) {
              match_gene[j] <- "Ambiguous"
            } else {
              match_gene[j] <- "Yes"
            }
          } else {
            match_gene[j] <- "No"
          }
          
        } # end loop over genes
        
        # combine information over genes
        if (!any(match_gene == "No")) {
          if (any(match_gene == "Ambiguous")) {
            match_vec[i] <- "Ambiguous"
          } else {
            match_vec[i] <- "Yes"
          }
        }
        
      } # end loop over data
      
      return(match_vec)
    }
  )
)


