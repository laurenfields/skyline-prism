# Skyline-PRISM (C#) dotnet-vNEXT Release Notes

Working draft for the next C# (.NET) release. Append entries as they land on the development branch;
rename to `RELEASE_NOTES_dotnet-v{version}.md` at release time - the release workflow publishes this file
as the GitHub Release description and fails if it is missing.

## New Features

- **`prism differential --report`: the quant report, headless.** The same report the Differential pane's
  **Quant report...** button writes (shipped in dotnet-v26.26.0), from the command's own flags, so it
  needs neither Skyline nor Windows:
  `prism differential -d output/ -g condition -a Control -b Disease --report --markers "EV markers (core)"`.
  `--markers` takes panels from the same set the Markers pane offers (saved lists plus the shipped
  panels), `--markers-group-by` picks their grouping column, and `--no-enrichment` skips g:Profiler on a
  machine with no internet access. The button and the command run one implementation, so a report
  clicked and a report typed agree: on the committed fixture, the report's `quant/differential.csv` is
  byte-identical to the results file the same command writes. A mistyped panel or column is refused
  before anything is written, and `--markers`, `--markers-group-by` or `--no-enrichment` without
  `--report` is refused rather than ignored.

## Bug Fixes

## Performance

## Breaking Changes
