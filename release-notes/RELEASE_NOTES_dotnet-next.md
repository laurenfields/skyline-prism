# Skyline-PRISM (C#) dotnet-vNEXT Release Notes

Working draft for the next C# (.NET) release. Append entries as they land on the development branch;
rename to `RELEASE_NOTES_dotnet-v{version}.md` at release time - the release workflow publishes this file
as the GitHub Release description and fails if it is missing.

## New Features

## Bug Fixes

- **A column of whole numbers is no longer fitted as a number just because it parses as one.** An
  integer-coded subject ID or batch used as a covariate was entered as one centered column - a
  straight-line effect of an arbitrary label - with nothing said. Now a column of numbers is
  categorical when its name contains a whole word such as patient, subject, id, batch, plate, cycle,
  set or run, or when it holds whole numbers with at most 10 distinct values; numeric otherwise.
  The type is shown beside each Adjust-for column (click to switch it) and can be set with
  `--covariate-type COLUMN=numeric|categorical`; results, `quant_parameters` (`covariate_types`) and
  recorded commands name the type each covariate was fitted as. A categorical covariate nested in
  the groups - a patient ID under a sex contrast - is dropped by name instead of failing with a
  generic rank-deficiency error, and a group-labelling column forced numeric is fitted with a
  warning. On a 63-sample ALS cohort, "adjusting" a sex contrast for patient IDs 1-15 had been
  reporting 78 peptides as a numeric covariate; it now reports that the covariate is nested and
  dropped, and the 86 unadjusted hits.

## Performance

## Breaking Changes

- **Covariate types changed for integer-coded columns.** A covariate column of whole numbers with
  at most 10 distinct values, or whose name contains patient, subject, donor, id, batch, plate,
  cycle, set or run, is now categorical by default where it used to be numeric, so an adjusted
  contrast using one gives different numbers than before. `--covariate-type COLUMN=numeric` (or
  clicking the type in the pane) restores the old model for that column.
