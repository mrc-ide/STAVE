# STAVE 2.0.4

* Updated the dependency on *variantstring* to v1.8.7. This version fixes an
  issue where `check_variant_string()` could fail with an empty error message,
  meaning that invalid variant strings in the counts table (or an invalid
  `target_variant` in `get_prevalence()`) now report an informative reason for
  the failure. It also includes some other bug fixes and efficiency
  improvements.
* Fixed a bug in the validation of the surveys table, where the compulsory
  `collection_day` field was not being checked correctly for `Date` class.
* Added the `sample_source` field to the surveys table, recording the sampling
  context from which samples were obtained (e.g. `clinical_passive`,
  `community_household`, `cohort`).
