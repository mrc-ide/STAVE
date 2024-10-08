% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{disallowed_variant_strings}
\alias{disallowed_variant_strings}
\title{Examples of invalid variant strings}
\usage{
disallowed_variant_strings()
}
\description{
Returns a data.frame of invalid variant strings, along with the reason for
being invalid. Can be useful in understanding what is the correct format.
Also, this data.frame is used in internal package testing, i.e. every row of
this data.frame should return an error.
}
