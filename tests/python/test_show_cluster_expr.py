"""Tests for the TO-GCN expression heatmap command-line program."""

from __future__ import annotations

import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

import numpy as np
import pandas as pd


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPOSITORY_ROOT / "scripts" / "python" / "show_cluster_expr.py"
SPEC = importlib.util.spec_from_file_location("show_cluster_expr", SCRIPT_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"Could not load module specification for {SCRIPT_PATH}")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CommandLineTests(unittest.TestCase):
    def test_parse_args_preserves_four_positionals_and_defaults_offset(self) -> None:
        args = MODULE.parse_args(
            ["expr.tsv", "levels.csv", "plot.png", "scaled.csv"]
        )

        self.assertEqual(args.expression_tsv, Path("expr.tsv"))
        self.assertEqual(args.level_csv, Path("levels.csv"))
        self.assertEqual(args.output_png, Path("plot.png"))
        self.assertEqual(args.output_csv, Path("scaled.csv"))
        self.assertEqual(args.level_offset, 1)

    def test_display_level_subtracts_configured_offset(self) -> None:
        self.assertEqual(MODULE.display_level(2, 1), "L1")
        self.assertEqual(MODULE.display_level(14, 1), "L13")
        self.assertEqual(MODULE.display_level(2, 0), "L2")


class InputValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temporary_directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary_directory.cleanup)
        self.directory = Path(self.temporary_directory.name)

    def write_inputs(
        self,
        expression_text: str,
        level_text: str,
    ) -> tuple[Path, Path]:
        expression_path = self.directory / "expression.tsv"
        level_path = self.directory / "levels.csv"
        expression_path.write_text(expression_text, encoding="utf-8")
        level_path.write_text(level_text, encoding="utf-8")
        return expression_path, level_path

    def test_load_inputs_returns_numeric_common_gene_tables(self) -> None:
        expression_path, level_path = self.write_inputs(
            "Gene ID\tSAM\tP1_P2\n"
            "gene_b\t2\t3\n"
            "gene_a\t1\t4\n"
            "expression_only\t5\t6\n",
            "TF gene ID,level in GCN\n"
            "gene_a,2\n"
            "gene_b,3\n"
            "level_only,4\n",
        )

        expression, levels = MODULE.load_inputs(expression_path, level_path)

        self.assertEqual(list(expression.index), ["gene_b", "gene_a"])
        self.assertEqual(list(expression.columns), ["SAM", "P1_P2"])
        self.assertEqual(expression.loc["gene_a", "P1_P2"], 4.0)
        self.assertEqual(list(levels.index), ["gene_a", "gene_b"])
        self.assertEqual(levels.loc["gene_a", "source_level"], 2)

    def test_load_inputs_rejects_duplicate_gene_ids(self) -> None:
        expression_path, level_path = self.write_inputs(
            "Gene ID\tSAM\n"
            "gene_a\t1\n"
            "gene_a\t2\n",
            "TF gene ID,level in GCN\n"
            "gene_a,2\n",
        )

        with self.assertRaisesRegex(ValueError, "duplicate gene IDs"):
            MODULE.load_inputs(expression_path, level_path)

    def test_load_inputs_rejects_nonnumeric_expression(self) -> None:
        expression_path, level_path = self.write_inputs(
            "Gene ID\tSAM\n"
            "gene_a\tnot-a-number\n",
            "TF gene ID,level in GCN\n"
            "gene_a,2\n",
        )

        with self.assertRaisesRegex(ValueError, "nonnumeric expression"):
            MODULE.load_inputs(expression_path, level_path)

    def test_load_inputs_rejects_nonnumeric_level(self) -> None:
        expression_path, level_path = self.write_inputs(
            "Gene ID\tSAM\n"
            "gene_a\t1\n",
            "TF gene ID,level in GCN\n"
            "gene_a,L2\n",
        )

        with self.assertRaisesRegex(ValueError, "numeric integer levels"):
            MODULE.load_inputs(expression_path, level_path)

    def test_load_inputs_rejects_empty_gene_intersection(self) -> None:
        expression_path, level_path = self.write_inputs(
            "Gene ID\tSAM\n"
            "expression_gene\t1\n",
            "TF gene ID,level in GCN\n"
            "level_gene,2\n",
        )

        with self.assertRaisesRegex(ValueError, "No gene IDs are shared"):
            MODULE.load_inputs(expression_path, level_path)


class ExpressionProcessingTests(unittest.TestCase):
    def test_order_genes_sorts_levels_and_handles_constant_and_single_gene(self) -> None:
        expression = pd.DataFrame(
            {
                "SAM": [3.0, 1.0, 5.0, 7.0],
                "P1_P2": [2.0, 2.0, 5.0, 8.0],
                "P3": [1.0, 3.0, 5.0, 9.0],
            },
            index=["gene_b", "gene_a", "constant_gene", "single_gene"],
        )
        levels = pd.DataFrame(
            {"source_level": [2, 2, 2, 3]},
            index=["gene_b", "gene_a", "constant_gene", "single_gene"],
        )

        ordered, groups = MODULE.order_genes_by_level(expression, levels)

        self.assertEqual([group[0] for group in groups], [2, 3])
        self.assertEqual(set(ordered[:2]), {"gene_a", "gene_b"})
        self.assertEqual(ordered[2:], ["constant_gene", "single_gene"])
        self.assertEqual(groups[0][1], ordered[:3])
        self.assertEqual(groups[1][1], ["single_gene"])

    def test_row_zscore_scales_variable_rows_and_zeroes_constant_rows(self) -> None:
        expression = pd.DataFrame(
            {"SAM": [1.0, 5.0], "P1_P2": [2.0, 5.0], "P3": [3.0, 5.0]},
            index=["variable_gene", "constant_gene"],
        )

        scaled = MODULE.row_zscore(expression)

        self.assertEqual(list(scaled.index), list(expression.index))
        self.assertEqual(list(scaled.columns), list(expression.columns))
        self.assertAlmostEqual(float(scaled.loc["variable_gene"].mean()), 0.0)
        self.assertAlmostEqual(float(scaled.loc["variable_gene"].std(ddof=0)), 1.0)
        np.testing.assert_array_equal(
            scaled.loc["constant_gene"].to_numpy(), np.zeros(3)
        )


class EndToEndTests(unittest.TestCase):
    def test_cli_creates_nonempty_heatmap_and_ordered_scaled_csv(self) -> None:
        with tempfile.TemporaryDirectory() as temporary_directory:
            directory = Path(temporary_directory)
            expression_path = directory / "expression.tsv"
            level_path = directory / "levels.csv"
            output_png = directory / "nested" / "expression.png"
            output_csv = directory / "nested" / "scaled.csv"
            expression_path.write_text(
                "Gene ID\tSAM\tP1_P2\tP3\n"
                "gene_a\t1\t2\t3\n"
                "gene_b\t3\t2\t1\n"
                "constant_gene\t5\t5\t5\n"
                "gene_c\t2\t3\t4\n"
                "expression_only\t8\t9\t10\n",
                encoding="utf-8",
            )
            level_path.write_text(
                "TF gene ID,level in GCN\n"
                "gene_a,2\n"
                "gene_b,2\n"
                "constant_gene,2\n"
                "gene_c,3\n"
                "level_only,4\n",
                encoding="utf-8",
            )

            completed = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT_PATH),
                    str(expression_path),
                    str(level_path),
                    str(output_png),
                    str(output_csv),
                ],
                capture_output=True,
                text=True,
                check=False,
            )

            self.assertEqual(completed.returncode, 0, completed.stderr)
            self.assertTrue(output_png.is_file())
            self.assertGreater(output_png.stat().st_size, 0)
            self.assertTrue(output_csv.is_file())
            self.assertNotIn(b"\r\n", output_csv.read_bytes())
            scaled = pd.read_csv(output_csv, index_col=0)
            self.assertEqual(
                set(scaled.index),
                {"gene_a", "gene_b", "constant_gene", "gene_c"},
            )
            self.assertEqual(list(scaled.columns), ["SAM", "P1_P2", "P3"])
            self.assertTrue(np.isfinite(scaled.to_numpy()).all())
            self.assertIn("Matched genes: 4", completed.stdout)
            self.assertIn("Expression-only genes excluded: 1", completed.stdout)
            self.assertIn("Level-only genes excluded: 1", completed.stdout)


if __name__ == "__main__":
    unittest.main()
