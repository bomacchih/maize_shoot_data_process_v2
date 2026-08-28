import importlib.util
import sys
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = (
    PROJECT_ROOT
    / "scripts"
    / "python"
    / "compare_supplementary_table9_markers.py"
)
SPEC = importlib.util.spec_from_file_location("compare_table9_markers", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
COMPARISON = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = COMPARISON
SPEC.loader.exec_module(COMPARISON)


def source_row(row_number, gene_id, cell_type):
    return {
        "source_row": row_number,
        "gene_id": gene_id,
        "cell_type": cell_type,
        "source_cell_type": cell_type.replace("_", " "),
        "celltype_id": "ontology-id",
    }


class SupplementaryTable9ComparisonTests(unittest.TestCase):
    def setUp(self):
        self.table_rows = [
            source_row(4, "Zm00001eb000001", "Type_One"),
            source_row(5, "Zm00001eb000001", "Type_One"),
            source_row(6, "Zm00001eb000002", "Type_One"),
            source_row(7, "Zm00001eb000003", "Type_One"),
            source_row(8, "Zm00001eb000003", "Type_Two"),
            source_row(9, "Zm00001eb000004", "Type_Two"),
        ]
        self.marker_rows = [
            {
                "gene_id": "Zm00001eb000001",
                "cell_type": "Type_One",
                "marker_rank": "1",
                "source_cell_type": "Type One",
            },
            {
                "gene_id": "Zm00001eb000002",
                "cell_type": "Type_One",
                "marker_rank": "2",
                "source_cell_type": "Type One",
            },
            {
                "gene_id": "Zm00001eb000004",
                "cell_type": "Type_Two",
                "marker_rank": "1",
                "source_cell_type": "Type Two",
            },
        ]

    def test_curated_exclusive_subset_passes(self):
        result = COMPARISON.compare(self.table_rows, self.marker_rows)
        self.assertTrue(result["comparison_passed"])
        self.assertEqual(result["table9_duplicate_pair_rows"], 1)
        self.assertEqual(result["table9_genes_shared_across_cell_types"], 1)
        self.assertEqual(result["marker_pairs_absent_from_table9"], 0)

    def test_shared_marker_fails(self):
        marker_rows = self.marker_rows + [
            {
                "gene_id": "Zm00001eb000003",
                "cell_type": "Type_One",
                "marker_rank": "3",
                "source_cell_type": "Type One",
            }
        ]
        result = COMPARISON.compare(self.table_rows, marker_rows)
        self.assertFalse(result["comparison_passed"])
        self.assertEqual(result["marker_pairs_not_celltype_exclusive_in_table9"], 1)


if __name__ == "__main__":
    unittest.main()
