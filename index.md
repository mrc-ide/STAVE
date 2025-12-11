
<img src="man/figures/stave.png" alt="STAVE Logo" style="max-height: 250px; margin-top:50px; margin-bottom:20px;">

STAVE stands for **S**patial-**T**emporal **A**ggregated **V**ariant **E**ncoding.

It is a lightweight data schema for organising aggregate counts of non-synonymous mutations, which is a data type commonly found in malaria molecular surveillance (MMS) studies. STAVE provides a consistent way to link genetic measurements to precise space–time coordinates, and to estimate the prevalence of particular genetic variants even when mixtures or ambiguous calls are present.

## What STAVE is designed to do

STAVE focuses on a small number of key tasks that commonly cause errors or inconsistencies in MMS workflows:

- **Provide a simple, relational data structure** for storing study metadata, survey-level space–time information, and aggregated genetic counts.
- **Encode haplotypes of non-synonymous mutations (i.e. changes in amino acid) in a compact, human-readable format** using the [variantstring](https://github.com/mrc-ide/variantstring) format.
- **Support robust prevalence estimation**, even when variants represent subsets of longer haplotypes or when mixed calls introduce uncertainty.

STAVE is deliberately minimal: it does not attempt to handle individual-level data, nucleotide-level data, or downstream spatial modelling. Instead, it aims to slot cleanly into broader analysis workflows by providing a reliable and unambiguous foundation on which those steps can build.

## Getting started

To learn how the STAVE data structure works, see [How it works](articles/core_design_principles.html).

When you’re ready to start working with the package, visit:

- [Installation](articles/installation.html) – how to install the package
- [Tutorials](articles/reading_in_data.html) – practical examples of common tasks

