
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
# wrap a character string at a fixed number of characters
#' @noRd
wrap_message <- function(x, width = 80) {
  paste(strwrap(x, width = width), collapse = "\n")
}
