
Files test the following qualities:

invalid_structure (expect fail):

invalid_structure_01: study_ids are not unique
invalid_structure_02: study_ids are not referenced in study_keys
invalid_structure_03: survey_ids are not unique within a study
invalid_structure_04: study_keys cannot be found in study_ids
invalid_structure_05: survey_ids are not referenced in survey_keys
invalid_structure_06: survey_keys cannot be found in survey_ids
invalid_structure_07: variant_strings are duplicated within a survey (although listed in different order)
invalid_structure_08: variant_num exceeds total_num
invalid_structure_09: total_num is not identical after grouping by survey-gene-locus
invalid_structure_10: variant_num sums to more than total_num within a survey-gene-locus combination

correct_structure (expect pass):

1. Tests that NAs in optional columns are treated correctly