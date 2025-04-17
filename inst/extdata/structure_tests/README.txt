
Files test the following qualities:

invalid_structure (expect fail):

invalid_structure_01: study_ids (in studies table) are not unique
invalid_structure_02: study_ids (in studies table) are not referenced in study_keys (in surveys table)
invalid_structure_03: combination of study_keys and survey_ids (in surveys table) are not unique
invalid_structure_04: study_keys (in surveys table) cannot be found in study_ids (in studies table)
invalid_structure_05: combination of study_keys and survey_ids (in surveys table) are not referenced in counts table
invalid_structure_06: combination of study_keys and survey_keys (in counts table) cannot be found in surveys table
invalid_structure_07: variant_strings are duplicated within a study-survey combination (although listed in different order)
invalid_structure_08: variant_num exceeds total_num
invalid_structure_09: total_num is not identical after grouping by study_id, survey_id, and gene-locus combination
invalid_structure_10: variant_num sums to more than total_num

correct_structure (expect pass):

correct_structure_01: checks that missing values are treated correctly