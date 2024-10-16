
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STAVE

[![master
checks](https://github.com/mrc-ide/STAVE/workflows/checks_main/badge.svg)](https://github.com/mrc-ide/STAVE/actions)
[![develop
checks](https://github.com/mrc-ide/STAVE/workflows/checks_develop/badge.svg)](https://github.com/mrc-ide/STAVE/actions)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/STAVE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/STAVE?branch=main)

STAVE is an R package for storing spatial-temporal genetic data. Data
are stored at the aggregate level (numerator and denominator) rather
than at the individual level. STAVE centres around one main class; data
are read into this class, and all data manipulation (e.g. calculating
prevalence) takes place using member functions. Although STAVE is
designed with *Plasmodium* genetic data in mind, in theory could be used
to store genetic information from other organisms.

## Why do we need this?

There are a few motivations for developing this package:

### 1. Efficient encoding of information

STAVE achieves a relatively efficient encoding through a relational
database structure. There are three tables that make up this database:

1.  **Studies**: These could be academic publications, reports, or other
    sources of information.
2.  **Surveys**: A survey is defined as a discrete instance of data
    collection, which includes information on geography (latitude and
    longitude) and collection time (day). There may be multiple surveys
    within a given study, for example different sampling sites or
    different collection periods within a single academic publication.
3.  **Counts**: This is where the genetic data are stored. Data include
    the name of the variant that we observed and how many times we
    observed it relative to a denominator count.

Separating the data into tables avoids duplication of information. For
example, the same author list applies to all surveys within a single
study, and therefore we do not want to repeat this information over all
surveys. Instead, we store this information at the study level but link
to the survey-level through a relational key, meaning we only need to
store the information once. This is both efficient in terms of memory
and reduces the opportunity for data entry mistakes.

### 2. Nailing down the input format

STAVE performs close to 100 different checks on the data when it is read
in. This includes everything from character formatting to checking that
prevalence cannot exceed 100%. Having this extremely pedantic series of
checks ensures that we know exactly what we are working with once it is
loaded, which makes life much easier when it comes to downstream
analysis.

### 3. Encoding genetic variants and calculating prevalence

Encode *Plasmodium* genetic data is quite tricky, and how we calculate
prevalence from these data is even trickier. Sequence data can be at the
single-codon level or spanning multiple codons, there can be
heterozygous calls (more than one allele observed) at some codons but
not others, heterozygous calls can either be phased or unphased, and
there can be different patterns of missing data between individuals that
should be reflected in different denominator counts. In STAVE we specify
a representation of the variant data that is flexible enough to
accommodate all of these situations.

A downside of this flexible encoding is that it is no longer trivial to
calculate the prevalence. This is kind of unavoidable - flexibly
encoding the information and easily calculating the prevalence are
somewhat opposing objectives. The approach in STAVE is to focus on the
flexible encoding so that we get everything into a nice clean structure,
and then to use member functions to calculate and return the prevalence
from the encoded data. Sometimes this prevalence calculation can be
quite complicated, for example in the case of an unphased multi-locus
variant it may only be possible to give lower and upper bounds on the
prevalence rather than calculating an exact value.

## Installation

You can install the development version of STAVE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("mrc-ide/STAVE")
```

## Example workflow

Input data take the form of three tables. We will use the pre-loaded
example data here, if you are using your own data then you need to match
this input format:

``` r
library(STAVE)

data("example_input")
View(example_input$studies)
View(example_input$surveys)
View(example_input$counts)
```

Create a STAVE object and append the data:

``` r
# create new object
s <- STAVE_object$new()

# append data using a member function
s$append_data(studies_dataframe = example_input$studies,
              surveys_dataframe = example_input$surveys,
              counts_dataframe = example_input$counts)
#> data correctly appended
```

``` r

# check how many studies are now loaded
s
#> Studies: 7
#> Surveys: 24
```

Once data are loaded, we can always view the different tables using
*get* functions. However, we cannot alter the values directly.

``` r
s$get_studies() |> head()
#> # A tibble: 6 × 6
#>   study_ID            study_name       study_type authors publication_year url  
#>   <chr>               <chr>            <chr>      <chr>              <dbl> <chr>
#> 1 wwarn_10297_Nelson  wwarn_10297_Nel… peer_revi… Nelson              1000 http…
#> 2 wwarn_10814_Dama    wwarn_10814_Dama peer_revi… Dama                1000 http…
#> 3 wwarn_10992_Mallick wwarn_10992_Mal… peer_revi… Mallick             1000 http…
#> 4 wwarn_11208_Kunasol wwarn_11208_Kun… peer_revi… Kunasol             1000 http…
#> 5 wwarn_11435_Henry   wwarn_11435_Hen… peer_revi… Henry               1000 http…
#> 6 wwarn_11720_Ould    wwarn_11720_Ould peer_revi… Ould                1000 http…
```

``` r
s$get_surveys() |> head()
#> # A tibble: 6 × 11
#>   study_key           survey_ID country_name site_name   lat   lon spatial_notes
#>   <chr>               <chr>     <chr>        <chr>     <dbl> <dbl> <chr>        
#> 1 wwarn_10297_Nelson  wwarn_10… Thailand     Kanchana…  15.3 98.5  wwarn lat an…
#> 2 wwarn_10814_Dama    wwarn_10… Mali         Koulikoro  12.6 -8.14 wwarn lat an…
#> 3 wwarn_10992_Mallick wwarn_10… India        Chhattis…  20.1 80.8  wwarn lat an…
#> 4 wwarn_10992_Mallick wwarn_10… India        Goa        15.3 74.1  wwarn lat an…
#> 5 wwarn_10992_Mallick wwarn_10… India        Gujarat    23.0 73.6  wwarn lat an…
#> 6 wwarn_10992_Mallick wwarn_10… India        Heilongj…  88.3 27.2  wwarn lat an…
#> # ℹ 4 more variables: collection_start <chr>, collection_end <chr>,
#> #   collection_day <chr>, time_notes <chr>
```

``` r
s$get_counts() |> head()
#> # A tibble: 6 × 4
#>   survey_key                           variant_string variant_num total_num
#>   <chr>                                <chr>                <dbl>     <dbl>
#> 1 wwarn_10297_Nelson_Sangkhlaburi_2002 mdr1:184:F              27        49
#> 2 wwarn_10297_Nelson_Sangkhlaburi_2002 mdr1:86:Y                4        49
#> 3 wwarn_10814_Dama_Bamako_2014         crt:76:T               130       170
#> 4 wwarn_10814_Dama_Bamako_2014         mdr1:86:Y               46       158
#> 5 wwarn_10992_Mallick_Assam_2002       crt:76:T                26        26
#> 6 wwarn_10992_Mallick_Assam_2002       mdr1:86:Y               19        25
```

We can calculate the prevalence of any variant using `get_prevalence()`:

``` r
p <- s$get_prevalence("mdr1:184:F")
View(p)
```

We can also return a list of all variants in the data:

``` r
s$get_variants()
#> [1] "crt:76:T"   "k13:469:F"  "k13:469:Y"  "k13:539:T"  "k13:580:Y" 
#> [6] "k13:675:V"  "mdr1:184:F" "mdr1:86:Y"
```

Finally, we can selectively drop studies from the data using their study
ID:

``` r
s$drop_study(drop_study_ID = "wwarn_10297_Nelson")
s
#> Studies: 6
#> Surveys: 23
```
