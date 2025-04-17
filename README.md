
<!-- README.md is generated from README.Rmd. Please edit that file -->

# STAVE

[![master
checks](https://github.com/mrc-ide/STAVE/workflows/checks_main/badge.svg)](https://github.com/mrc-ide/STAVE/actions)
[![develop
checks](https://github.com/mrc-ide/STAVE/workflows/checks_develop/badge.svg)](https://github.com/mrc-ide/STAVE/actions)
[![Codecov test
coverage](https://codecov.io/gh/mrc-ide/STAVE/branch/main/graph/badge.svg)](https://app.codecov.io/gh/mrc-ide/STAVE?branch=main)

STAVE is an R package with three basic functionalities:

1.  **Stores aggregate genetic data (i.e. numerator and denominator) in
    [variant string format](https://github.com/mrc-ide/variantstring)**.
    This is a convenient format when dealing with amino acid level data,
    as is common for drug resistance markers. It can handle single
    codons as well as multi-locus haplotypes.
2.  **Links genetic data to a specific point in space and time**. This
    avoids any issues with country or admin names changing or being
    spelled differently.
3.  **Calculates prevalence from the encoded data**. This is not as
    simple as it sounds, for example the variant of interest may be a
    subset of an encoded haplotype, or there may be mixed (heterozygous)
    calls in the data.

All documentation, including installation instructions and tutorials,
can be found on [the STAVE website](https://mrc-ide.github.io/STAVE/).

# Version History

The current version of the software is 1.1.0, released 17 Apr 2025.
