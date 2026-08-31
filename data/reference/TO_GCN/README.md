# TO-GCN reference inputs

This directory stores the small reference tables required to reconstruct the
time-ordered gene coexpression network (TO-GCN). Large Seurat objects remain in
`data/processed/` and are distributed separately from GitHub.

Provide one of the following transcription-factor annotations:

- `maize_TF_annotation.csv` with columns `gene_id`, `gene_name`, and
  `tf_family`; the provided table also has the optional `archived_level`
  comparison column; or
- `gene_expressions_v3.xlsx`, the recovered analysis workbook, whose
  `TF_genes` sheet contains the maize TF annotation and archived L1-L13 level
  assignments.

The archived level assignments are used only as a clearly labeled fallback for
reconstructing the published demonstration figure when fresh `Node_level.csv`
output is absent. New datasets must be processed with the TO-GCN executables;
their levels must not be copied from the archived workbook.

Optional local GO enrichment input:

- `maize_gene_to_GO.csv` with columns `gene_id` and `GO_ID`. Multiple GO IDs in
  one record may be separated by spaces, commas, semicolons, or vertical bars.

Without the optional GO mapping, the postprocessing script still exports one
AgriGO-compatible gene list for every TO-GCN level and an expressed-gene
background list.

The original single-time-series software is obtained from
<https://github.com/petitmingchang/TO-GCN/tree/master/Single_Time-series_data>.
The compilation script records the exact source commit used for each run.
