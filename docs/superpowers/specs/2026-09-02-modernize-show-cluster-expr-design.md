# Modernize `show_cluster_expr.py` Design

## Objective

Replace the legacy positional-script implementation with a reliable, documented,
testable command-line program while preserving its four required positional
arguments:

```text
python show_cluster_expr.py EXPRESSION_TSV TF_LEVEL_CSV OUTPUT_PNG OUTPUT_CSV
```

The repository copy will live at `scripts/python/show_cluster_expr.py`. The
supplied raw inputs will live under `results/tables/11_TO_GCN/input/`.

## Command-line interface

The program will use `argparse` and accept these required positional paths:

1. expression matrix: tab-separated, with gene IDs in the first column;
2. TF-level table: comma-separated, with TF gene IDs in the first column and
   numeric source levels in the last column;
3. output PNG path; and
4. output rescaled-expression CSV path.

It will also accept `--level-offset INTEGER`, defaulting to `1`. Displayed
levels are calculated as `source level - level offset`. With the supplied
`TF_level.csv`, source levels 2-14 are therefore displayed as L1-L13. Setting
`--level-offset 0` preserves the source labels.

## Processing design

The script will be organized into focused functions for argument parsing, input
loading and validation, within-level gene ordering, row-wise scaling, plotting,
and orchestration through `main()`.

Input validation will reject missing files, malformed tables, duplicate gene
IDs, nonnumeric expression values or levels, an empty gene intersection, and
non-finite values. Level groups will follow ascending numeric source-level
order. Genes absent from either input will be reported and excluded without
changing the input files.

Within each level, genes will be ordered by average-linkage hierarchical
clustering of Pearson-correlation distance. The implementation will pass the
condensed distance vector expected by SciPy rather than treating a square
distance matrix as observations. A one-gene group requires no clustering.
Constant-expression genes have undefined correlation distance, so they will be
ordered deterministically after variable genes and will receive zeroes during
row-wise z-score scaling.

The rescaled CSV will contain genes in exactly the plotted order and retain the
five expression-column names. Output parent directories will be created when
needed. The plot will use a noninteractive Matplotlib backend, preserve the
current heatmap dimensions, color range, and high-resolution PNG output, draw
boundaries between levels, and label groups using the offset-adjusted levels.
Fragile position-based annotations from the legacy script will be removed.

## Repository documentation and artifacts

`docs/11_TO_GCN_embryonic_leaf.md` will gain a dedicated heatmap step explaining
the input formats, default level offset, exact repository-root command, outputs,
and interpretation of z-scores. The root `README.md` step-11 entry will point to
the Python heatmap script as part of the TO-GCN implementation.

The supplied `tf_genes.txt` and `TF_level.csv` will be copied unchanged to
`results/tables/11_TO_GCN/input/`. The verified command will write
`results/figures/11_TO_GCN/expression.png` and
`results/tables/11_TO_GCN/postprocess/rescaled_expr.csv`. These derived files
will be committed so the documented example has inspectable outputs.

## Testing and verification

Tests in `tests/python/test_show_cluster_expr.py` will be written before the new
implementation. They will cover argument parsing and the default offset,
correct level-label transformation, deterministic handling of single and
constant genes, input validation, ordered scaled-output content, and an
end-to-end CLI run that creates both outputs.

Verification will include the focused test file, the repository's existing
Python tests, `--help`, the exact full-data command, structural inspection of the
resulting CSV, and visual inspection of `expression.png`.

## Reviewer-response accuracy

The revised response will explicitly flag two unresolved protocol facts for
coauthor action instead of presenting them as settled:

- the verified Visium workflow applies curated metadata membership but no
  numeric spot-level cutoffs for detected genes, total UMIs, mitochondrial
  percentage, or chloroplast percentage; and
- the supplied `TF_level.csv` contains 1,264 TF assignments across 13 source
  levels (2-14), displayed as L1-L13 by the plotting script's default offset.

All other reported parameters will be taken directly from the deposited scripts
and records: 30 PCs, seed 1234, Harmony dimensions 1-30, resolution 2, UMAP
neighbors 30/minimum distance 0.3/cosine metric, and Wilcoxon positive-marker
testing with `min.pct = 0.25`, `logfc.threshold = 0.25`, and adjusted P < 0.05.
