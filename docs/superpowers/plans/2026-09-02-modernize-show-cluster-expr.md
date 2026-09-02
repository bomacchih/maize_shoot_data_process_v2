# TO-GCN Heatmap Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Deliver a tested `argparse`-based TO-GCN TF heatmap script, bundled example inputs and outputs, updated documentation, and an evidence-based reviewer response.

**Architecture:** A single importable Python module separates CLI parsing, validation/loading, level-wise correlation clustering, z-score scaling, plotting, and orchestration. Standard `unittest` tests import the module directly and exercise pure functions plus one subprocess-level CLI run.

**Tech Stack:** Python 3, argparse, pathlib, NumPy, pandas, SciPy, Matplotlib, unittest, Git.

---

### Task 1: Specify and implement the CLI and input contract

**Files:**
- Create: `tests/python/test_show_cluster_expr.py`
- Create: `scripts/python/show_cluster_expr.py`

- [ ] **Step 1: Write the failing CLI/default test**

Create a unittest module that loads `scripts/python/show_cluster_expr.py` with
`importlib.util`, then asserts:

```python
args = MODULE.parse_args(["expr.tsv", "levels.csv", "plot.png", "scaled.csv"])
self.assertEqual(args.level_offset, 1)
self.assertEqual(args.expression_tsv, Path("expr.tsv"))
```

Also assert that `display_level(2, 1) == "L1"`,
`display_level(14, 1) == "L13"`, and `display_level(2, 0) == "L2"`.

- [ ] **Step 2: Run the focused test to verify RED**

Run:

```powershell
python -m unittest tests.python.test_show_cluster_expr -v
```

Expected: import failure because `scripts/python/show_cluster_expr.py` does not
exist.

- [ ] **Step 3: Implement the CLI skeleton**

Create the script with an `ArgumentParser` containing four `Path` positional
arguments and:

```python
parser.add_argument(
    "--level-offset",
    type=int,
    default=1,
    help="Integer subtracted from source levels when labels are displayed "
         "(default: 1, so source levels 2-14 become L1-L13).",
)
```

Implement `display_level(source_level: int, level_offset: int) -> str` as
`f"L{source_level - level_offset}"`, and define `main(argv=None) -> int` with a
standard `if __name__ == "__main__"` exit guard.

- [ ] **Step 4: Verify GREEN**

Run the focused test and confirm both CLI/default and offset assertions pass.

### Task 2: Implement reliable loading, ordering, scaling, and plotting through TDD

**Files:**
- Modify: `tests/python/test_show_cluster_expr.py`
- Modify: `scripts/python/show_cluster_expr.py`

- [ ] **Step 1: Add failing validation tests**

Use `tempfile.TemporaryDirectory` to create small TSV/CSV fixtures. Test that
`load_inputs()` returns numeric expression data and numeric levels for valid
files, and raises `ValueError` containing stable messages for duplicate gene
IDs, nonnumeric values, and an empty gene intersection.

- [ ] **Step 2: Verify validation tests fail for missing implementation**

Run the focused unittest module and confirm failures reference `load_inputs`.

- [ ] **Step 3: Implement validation/loading**

Implement:

```python
def load_inputs(expression_path: Path, levels_path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    expression = pd.read_csv(expression_path, sep="\t", index_col=0)
    levels = pd.read_csv(levels_path)
```

Require at least one expression column and two level-table columns, nonblank
unique gene IDs, finite numeric expression values, integer-valued numeric source
levels, unique TF IDs, and a nonempty intersection. Preserve the level table's
column names while normalizing its first and last columns internally.

- [ ] **Step 4: Add failing ordering/scaling tests**

Test that `order_genes_by_level()` sorts numeric levels, retains all common
genes once, clusters variable genes deterministically, places constant genes
after variable genes within a level, and handles a single-gene level. Test that
`row_zscore()` produces per-row mean zero and converts a constant row to zeros.

- [ ] **Step 5: Implement ordering and scaling**

For each source level, split constant and variable rows. For two or more variable
genes, calculate:

```python
distances = scipy.spatial.distance.pdist(values, metric="correlation")
tree = scipy.cluster.hierarchy.linkage(distances, method="average", optimal_ordering=True)
ordered = [gene_ids[i] for i in scipy.cluster.hierarchy.leaves_list(tree)]
```

Append sorted constant gene IDs, and concatenate levels in numeric order.
Implement row z-scores explicitly with NumPy so zero-standard-deviation rows
become zero without warnings.

- [ ] **Step 6: Add a failing end-to-end CLI test**

Run the script with a fixture containing source levels 2 and 3. Assert exit code
zero, nonempty PNG and CSV files, the expected gene order/index, finite scaled
values, and output-column equality with the input expression columns.

- [ ] **Step 7: Implement plotting and orchestration**

Use Matplotlib's `Agg` backend and `imshow` for a 16 x 4.5 inch heatmap of the
transposed scaled matrix with `vmin=-1.5`, `vmax=1.5`, and `seismic`. Draw level
boundaries, center offset-adjusted labels under each group, label the color bar
`Z-score`, create output parents, write the ordered scaled CSV, and print a
concise completion summary including matched and excluded gene counts.

- [ ] **Step 8: Run focused tests to GREEN**

Run:

```powershell
python -m unittest tests.python.test_show_cluster_expr -v
```

Expected: all heatmap tests pass with no runtime warnings.

### Task 3: Add supplied data, documentation, and verified full-data artifacts

**Files:**
- Create: `results/tables/11_TO_GCN/input/tf_genes.txt`
- Create: `results/tables/11_TO_GCN/input/TF_level.csv`
- Create: `results/tables/11_TO_GCN/postprocess/rescaled_expr.csv`
- Create: `results/figures/11_TO_GCN/expression.png`
- Modify: `docs/11_TO_GCN_embryonic_leaf.md`
- Modify: `README.md`

- [ ] **Step 1: Copy raw inputs unchanged**

Copy the workspace files byte-for-byte and verify SHA-256 equality between each
source and repository destination.

- [ ] **Step 2: Update the TO-GCN workflow documentation**

Add a heatmap subsection with the exact repository-root command:

```powershell
python scripts/python/show_cluster_expr.py `
  results/tables/11_TO_GCN/input/tf_genes.txt `
  results/tables/11_TO_GCN/input/TF_level.csv `
  results/figures/11_TO_GCN/expression.png `
  results/tables/11_TO_GCN/postprocess/rescaled_expr.csv
```

State that `--level-offset` defaults to 1, source levels 2-14 display as L1-L13,
and `--level-offset 0` preserves source numbering. Document the actual 1,264-row
intersection and explain that scaling is gene-wise across the five domains.

- [ ] **Step 3: Update root README step 11**

Add `scripts/python/show_cluster_expr.py` to the step-11 implementation links
and mention the complete command in the TO-GCN report.

- [ ] **Step 4: Run the full-data command**

Run the documented command. Expected outputs are a nonempty PNG and a 1,264-row
CSV with five numeric columns and no non-finite values.

- [ ] **Step 5: Inspect the generated figure**

Open `results/figures/11_TO_GCN/expression.png` and verify readable L1-L13
labels, correct five-domain y-axis labels, clean boundaries, and no clipping.

### Task 4: Repository-wide verification, response drafting, commit, and push

**Files:**
- Modify as required only if verification identifies an in-scope defect.

- [ ] **Step 1: Run all Python unit tests**

Run:

```powershell
python -m unittest discover -s tests/python -p "test_*.py" -v
```

Expected: zero failures and zero errors.

- [ ] **Step 2: Verify CLI and repository integrity**

Run `python scripts/python/show_cluster_expr.py --help`, `git diff --check`,
input/output structural checks, and `git status --short`.

- [ ] **Step 3: Draft the reviewer response from verified evidence**

Report 30 PCs; seed 1234; Harmony dimensions 1-30 while highlighting that no
explicit Harmony grouping variable is supplied; resolution 2; UMAP 30
neighbors, minimum distance 0.3, cosine metric; and positive Wilcoxon markers
with `min.pct=0.25`, `logfc.threshold=0.25`, adjusted P < 0.05. Highlight for
coauthor action that numeric spot-level QC thresholds are absent. Correct the
TO-GCN input statement to 1,264 TFs across 13 source levels, displayed as L1-L13
through the default offset.

- [ ] **Step 4: Commit the implementation**

Stage only the plan, script, tests, supplied inputs, documented generated
outputs, README, and TO-GCN report. Commit with:

```text
feat: modernize TO-GCN expression heatmap
```

- [ ] **Step 5: Push and confirm remote state**

Push `main` to `origin`, then compare local `HEAD` with `origin/main`. Expected:
the commit hashes match and the worktree is clean.
