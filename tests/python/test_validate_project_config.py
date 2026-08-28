import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[2]
MODULE_PATH = PROJECT_ROOT / "scripts" / "python" / "validate_project_config.py"
SPEC = importlib.util.spec_from_file_location("validate_project_config", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
VALIDATOR = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = VALIDATOR
SPEC.loader.exec_module(VALIDATOR)


class ConfigurationValidatorTests(unittest.TestCase):
    def test_current_configuration_parses(self):
        config = VALIDATOR.load_yaml_subset(
            PROJECT_ROOT / "config" / "current_state.yaml"
        )
        self.assertEqual(config["schema_version"], 1)
        self.assertEqual(config["study"]["sample_ids"], VALIDATOR.EXPECTED_SAMPLES)
        self.assertEqual(config["developmental_trends_and_go"]["trend_clusters"], 7)

    def test_invalid_indentation_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.yaml"
            path.write_text("root:\n   key: 1\n", encoding="utf-8")
            with self.assertRaises(VALIDATOR.ValidationInputError):
                VALIDATOR.load_yaml_subset(path)

    def test_duplicate_key_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "bad.yaml"
            path.write_text("root:\n  key: 1\n  key: 2\n", encoding="utf-8")
            with self.assertRaises(VALIDATOR.ValidationInputError):
                VALIDATOR.load_yaml_subset(path)

    def test_schema_and_script_snapshot_have_no_failures(self):
        config = VALIDATOR.load_yaml_subset(
            PROJECT_ROOT / "config" / "current_state.yaml"
        )
        reporter = VALIDATOR.Reporter()
        VALIDATOR.validate_config_schema(config, reporter)
        VALIDATOR.validate_sample_manifest(PROJECT_ROOT, config, reporter)
        VALIDATOR.validate_script_consistency(PROJECT_ROOT, config, reporter)
        failures = [item for item in reporter.findings if item.severity == "FAIL"]
        self.assertEqual(failures, [])


if __name__ == "__main__":
    unittest.main()
