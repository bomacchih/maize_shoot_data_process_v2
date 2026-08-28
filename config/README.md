# Configuration snapshot

This directory records the analysis settings that are currently embedded in
the repository scripts. It separates portable scientific parameters from
machine-specific paths.

## Files

- `current_state.yaml` — consolidated snapshot of parameters, input status,
  and unresolved reproducibility gaps found in the current scripts.
- `paths.example.yaml` — template for external reference and input locations.
- `spaceranger_samples.csv` — sample order and the non-path Space Ranger
  settings recoverable from the current shell script.

## Important limitation

The analysis scripts do **not** currently read `current_state.yaml`. It is an
auditable description of the current implementation, not yet a control file.
Until the scripts are refactored to load this configuration, changing a value
here will not change an analysis run.

Fields marked `null`, `not_explicit`, or `missing_blocking` are intentionally
unresolved. They should not be replaced with guessed values. In particular,
the current Visium integration script does not define numeric spot-level QC
cutoffs, does not run clustering, and imports `harmony_clusters` from curated
metadata.

## Local paths

Copy `paths.example.yaml` to `local_paths.yaml` and fill in paths on the machine
that will run Space Ranger or rebuild raw-data workflows. `local_paths.yaml` is
excluded from Git because it may contain private filesystem locations.

Use project-relative paths whenever the input is stored inside this repository.
Do not commit credentials, private download links, reviewer links, or access
tokens.

## Published data

The public analysis files are deposited at
<https://zenodo.org/records/22058284> (DOI `10.5281/zenodo.22058284`). Download
the record files into `data/processed/` without renaming them. See
`data/processed/README.md` for the file groups and limitations.

The record contains 14 web summaries, 14 `.cloupe` files, 14 Velocyto `.loom`
files, the 14 per-sample Seurat objects, and the combined Visium and scRNA RDS
objects. It does not contain complete Space Ranger `outs` folders: the matrix
H5 and spatial image/position/scalefactor files needed to rebuild the individual
Seurat objects are not present in the record.

The combined `data/processed/sc_merged_filter_SCT2_inte.rds` is sufficient for
downstream scRNA annotation. The "raw scRNA matrices" are the three 10x files
(`matrix.mtx.gz`, `barcodes.tsv.gz`, and `features.tsv.gz`) for each of
SRR11943512 and SRR11943513; they are needed only to regenerate that combined
RDS from counts.

The curated SCINA marker source is `marker_list2.rds`, a named list derived
from Supplementary Table 9. Convert it with
`scripts/R/10_scRNA_reference_integration/00_prepare_SCINA_marker_table.R`.
The generated CSV retains the original cell-type label, uses a punctuation-safe
analysis label, and records each gene's position within its source list as
`marker_rank`. This preserves source order without claiming an effect-size
ranking that the RDS does not provide.

## Validation

From the repository root, run:

```bash
python scripts/python/validate_project_config.py
```

The validator uses only the Python standard library. It checks this snapshot,
the sample manifest, metadata integrity, files declared present, and important
script assignments. Add `--strict` when preparing a release so unresolved
warnings produce a nonzero exit status. See `tests/README.md` for the complete
test workflow, including Seurat object validation.
