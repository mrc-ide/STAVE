
#------------------------------------------------
# check that character string is a valid identifier. This means:
# -	It consists only of English letters (upper or lower case), numbers (0-9), or underscores (_).
# -	The first character cannot be a number or an underscore.
#' @noRd

assert_valid_string <- function(input_strings, message_name) {
  is_valid <- TRUE
  
  # Ensure the input is a character vector
  if (!is.character(input_strings)) {
    is_valid <- FALSE
  }
  
  # Iterate over each string in the vector and check each one
  for (input_string in input_strings) {
    
    # Check if the first character is a letter (upper or lower case)
    if (!grepl("^[A-Za-z]", input_string)) {
      is_valid <- FALSE
    }
    
    # Check if the entire string contains only allowed characters (letters, numbers, underscore)
    if (!grepl("^[A-Za-z0-9_]+$", input_string)) {
      is_valid <- FALSE
    }
  }
  
  # stop if invalid
  if (!is_valid) {
    stop(paste(c(sprintf("%s is not a valid identifier. It must:", message_name),
                 "- consist only of English letters (upper or lower case), numbers (0-9), or underscores (_).",
                 "- The first character cannot be a number or an underscore."), collapse = "\n"))
  }
}

#------------------------------------------------
# check that date_string conforms to YYYY-MM-DD format
#' @importFrom stringr str_detect
#' @importFrom lubridate ymd
#' @noRd
assert_valid_ymd <- function(date_string) {
  is_valid <- TRUE
  
  for (x in date_string) {
    
    # check string formatted correctly
    if (!stringr::str_detect(x, "^\\d{4}-\\d{1,2}-\\d{1,2}$")) {
      is_valid <- FALSE
    }
    
    # check can be read in as date
    if (is.na(lubridate::ymd(x, quiet = TRUE))) {
      is_valid <- FALSE
    }
    
  }
  
  # stop if invalid
  if (!is_valid) {
    stop("dates must be formatted in YYYY-MM-DD format")
  }
}

#------------------------------------------------
# check that date is greater than or equal to a value
#' @noRd
assert_date_greq <- function(x, y, message = NULL,
                             name_x = paste(deparse(substitute(x)), collapse = ""),
                             name_y = nice_format(y)) {
  
  # default message
  if (is.null(message)) {
    message <- sprintf("%s must be greater than or equal to %s", name_x, name_y)
  }
  
  assert_class(x, "Date")
  assert_class(y, "Date")
  assert_in(length(y), c(1, length(x)), message = message)
  if (!all(x >= y)) {
    stop(message, call. = FALSE)
  }
  return(TRUE)
}
#------------------------------------------------
# check that date is between two bounds
#' @noRd
assert_date_bounded <- function(d, left, right) {
  assert_class(d, "Date")
  assert_class(left, "Date")
  assert_class(right, "Date")
  
  if (any((d < left) | (d > right))) {
    stop("date must be within bounds")
  }
  
}

#------------------------------------------------
# drop amino acid sequence from a variant string
#' @noRd
drop_amino <- function(variant_string, sort_gene_name = FALSE) {
  
  # split by semicolon, drop everything from the last colon onwards, and paste
  # genes back together separated by semicolon
  ret <- mapply(function(x) {
    z <- sub(":[^:]*$", "", x)
    if (sort_gene_name) {
      z <- sort(z)
    }
    paste(z, collapse = ";")
  }, strsplit(variant_string, ";"))
  
  return(ret)
}

#------------------------------------------------
# reorder a variant string so genes are listed alphabetically
#' @noRd
sort_gene_name <- function(variant_string) {
  mapply(function(x) paste(sort(x), collapse = ";"), strsplit(variant_string, ";"))
}

#------------------------------------------------
# take a variant string and expand the elements into a nested list where they
# are easier to operate on.
# Can take inputs with amino acid sequence (e.g. "crt:72:C;dhps:540:E") or
# without (e.g. "crt:72;dhps:540")
#' @noRd
variant_to_list <- function(variant_string) {
  
  # this only works on single elements
  assert_single_string(variant_string)
  
  ret <- list()
  s1 <- strsplit(variant_string, split = ";")[[1]]
  for (i in seq_along(s1)) {
    s2 <- strsplit(s1[i], ":")[[1]]
    ret[[s2[1]]] <- list(pos = as.numeric(strsplit(s2[2], split = "_")[[1]]))
    if (length(s2) == 3) {
      ret[[s2[1]]]$amino <- unlist(strsplit(gsub("_", "", s2[3]), "(?<=[A-Z])(?=[^/|])", perl = TRUE))
    }
  }
  
  return(ret)
}

#------------------------------------------------
#' @title String of allowed amino acids
#'
#' @description 
#' Returns a data.frame of allowed amino acid single-letter codes. These come
#' from IUPAC (International Union of Pure and Applied Chemistry),
#' \href{https://www.bioinformatics.org/sms/iupac.html}{see here} for details.
#' Variants defined in the \code{variant_string} column of the counts data.frame
#' must be one of these single-letter codes, or otherwise the underscore (_),
#' pipe (|) or forward slash (/) symbols. Underscores will be removed but are
#' allowed for improved readability. The pipe and forward slash symbols
#' represent phased or unphased mixed calls, respectively.
#' 
#' @export

allowed_amino_acids <- function() {
  
  # get the full path to the file in the extdata directory
  file_path <- system.file("extdata", "IUPAC_codes.rds", package = "STAVE")
  
  # read the data and return
  ret <- readRDS(file_path)
  return(ret)
}

#------------------------------------------------
#' @title Examples of invalid variant strings
#'
#' @description 
#' Returns a data.frame of invalid variant strings, along with the reason for
#' being invalid. Can be useful in understanding what is the correct format.
#' Also, this data.frame is used in internal package testing, i.e. every row of
#' this data.frame should return an error.
#' 
#' @export

disallowed_variant_strings <- function() {
  
  # get the full path to the file in the extdata directory
  file_path <- system.file("extdata", "disallowed_variant_strings.rds", package = "STAVE")
  
  # read the data and return
  ret <- readRDS(file_path)
  return(ret)
}

#------------------------------------------------
# wrap a character string at a fixed number of characters
#' @noRd
wrap_message <- function(x, width = 80) {
  paste(strwrap(x, width = width), collapse = "\n")
}
