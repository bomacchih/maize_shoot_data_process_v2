import importlib.util
import sys
import unittest
from pathlib import Path


MODULE_PATH = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "python"
    / "validate_public_readiness.py"
)
SPEC = importlib.util.spec_from_file_location("validate_public_readiness", MODULE_PATH)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


class PublicReadinessTextScanTests(unittest.TestCase):
    def test_secret_location_is_reported_without_value(self):
        secret_value = "ghp_" + "A" * 36
        findings = MODULE.classify_text_findings(
            f"TOKEN={secret_value}\n", "example.txt"
        )
        self.assertEqual(findings["GITHUB_TOKEN"], {"example.txt:1"})
        self.assertNotIn(secret_value, repr(findings))

    def test_placeholder_assignment_is_not_a_secret(self):
        findings = MODULE.classify_text_findings(
            "api_key = your_api_key_placeholder\n", "example.yaml"
        )
        self.assertNotIn("GENERIC_SECRET_ASSIGNMENT", findings)

    def test_personal_home_path_is_reported(self):
        findings = MODULE.classify_text_findings(
            "input=/home/researcher/private/data\n", "workflow.sh"
        )
        self.assertEqual(
            findings["PERSONAL_ABSOLUTE_PATH"], {"workflow.sh:1"}
        )

    def test_generic_example_home_path_is_allowed(self):
        findings = MODULE.classify_text_findings(
            "input=/home/jdoe/example/data\n", "workflow.sh"
        )
        self.assertNotIn("PERSONAL_ABSOLUTE_PATH", findings)

    def test_r_s4_slot_is_not_reported_as_email(self):
        findings = MODULE.classify_text_findings(
            "is.na(object@meta.data[, columns, drop = FALSE])\n",
            "workflow.R",
        )
        self.assertNotIn("EMAIL_ADDRESS", findings)


class PublicReadinessApprovedDataTests(unittest.TestCase):
    def setUp(self):
        self.repository_root = Path(__file__).resolve().parents[2]
        self.relative_path = (
            "data/reference/scRNA_reference/marker_list2.rds"
        )
        self.content = (
            self.repository_root / self.relative_path
        ).read_bytes()

    def test_reviewed_marker_rds_is_approved(self):
        self.assertTrue(
            MODULE.is_approved_public_data_content(
                self.relative_path, self.content
            )
        )

    def test_changed_marker_rds_is_blocked(self):
        self.assertFalse(
            MODULE.is_approved_public_data_content(
                self.relative_path, self.content + b"changed"
            )
        )

    def test_same_content_at_another_path_is_blocked(self):
        self.assertFalse(
            MODULE.is_approved_public_data_content(
                "data/reference/scRNA_reference/unreviewed.rds",
                self.content,
            )
        )


if __name__ == "__main__":
    unittest.main()
