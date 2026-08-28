#!/usr/bin/env python3
"""Validate the maize shoot configuration, inputs, and script consistency.

This validator intentionally uses only the Python standard library. The project
configuration uses a conservative YAML subset (nested mappings and scalar
lists), which is parsed here without requiring PyYAML.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Iterable


EXPECTED_SAMPLES = [
    "UL01", "UL02", "UL04",
    "VR01", "VR02", "VR03", "VR04",
    "DQ01", "DQ02", "DQ03", "DQ04", "DQ06", "DQ07", "DQ08",
]
EXPECTED_DOMAINS = ["SAM", "P1_P2", "P3", "P4", "P5", "coleoptile", "co_v"]
REQUIRED_METADATA_COLUMNS = {
    "Barcode", "sample_id", "section_id", "domains", "harmony_clusters"
}


class ValidationInputError(ValueError):
    """Raised when a configuration file cannot be interpreted safely."""


@dataclass
class Finding:
    severity: str
    code: str
    message: str


class Reporter:
    def __init__(self) -> None:
        self.findings: list[Finding] = []

    def passed(self, code: str, message: str) -> None:
        self.findings.append(Finding("PASS", code, message))

    def warn(self, code: str, message: str) -> None:
        self.findings.append(Finding("WARN", code, message))

    def fail(self, code: str, message: str) -> None:
        self.findings.append(Finding("FAIL", code, message))

    def count(self, severity: str) -> int:
        return sum(item.severity == severity for item in self.findings)

    def print_report(self) -> None:
        for item in self.findings:
            print(f"[{item.severity}] {item.code}: {item.message}")
        print(
            "\nSummary: "
            f"{self.count('PASS')} passed, "
            f"{self.count('WARN')} warning(s), "
            f"{self.count('FAIL')} failure(s)."
        )


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value == "":
        return ""
    if value.startswith('"') and value.endswith('"'):
        try:
            return json.loads(value)
        except json.JSONDecodeError as error:
            raise ValidationInputError(f"Invalid quoted scalar: {value}") from error
    if value.startswith("'") and value.endswith("'"):
        return value[1:-1].replace("''", "'")
    lowered = value.lower()
    if lowered in {"null", "~"}:
        return None
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    if re.fullmatch(r"[-+]?(?:\d+\.\d*|\d*\.\d+)(?:[eE][-+]?\d+)?", value):
        return float(value)
    if re.fullmatch(r"[-+]?\d+[eE][-+]?\d+", value):
        return float(value)
    return value


def load_yaml_subset(path: Path) -> dict[str, Any]:
    """Load the mapping/list/scalar YAML subset used by current_state.yaml."""
    logical_lines: list[tuple[int, int, str]] = []
    for line_number, raw_line in enumerate(
        path.read_text(encoding="utf-8-sig").splitlines(), start=1
    ):
        if "\t" in raw_line:
            raise ValidationInputError(f"{path}:{line_number}: tabs are not allowed")
        stripped = raw_line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        indent = len(raw_line) - len(raw_line.lstrip(" "))
        if indent % 2:
            raise ValidationInputError(
                f"{path}:{line_number}: indentation must use multiples of two spaces"
            )
        logical_lines.append((line_number, indent, stripped))

    if not logical_lines:
        raise ValidationInputError(f"{path}: configuration is empty")

    def parse_block(position: int, indent: int) -> tuple[Any, int]:
        line_number, actual_indent, content = logical_lines[position]
        if actual_indent != indent:
            raise ValidationInputError(
                f"{path}:{line_number}: unexpected indentation {actual_indent}"
            )
        is_list = content.startswith("- ")
        container: Any = [] if is_list else {}

        while position < len(logical_lines):
            line_number, actual_indent, content = logical_lines[position]
            if actual_indent < indent:
                break
            if actual_indent > indent:
                raise ValidationInputError(
                    f"{path}:{line_number}: unexpected nested content"
                )

            if is_list:
                if not content.startswith("- "):
                    raise ValidationInputError(
                        f"{path}:{line_number}: cannot mix a list and mapping at one level"
                    )
                item = content[2:].strip()
                if not item:
                    raise ValidationInputError(
                        f"{path}:{line_number}: nested list items are not supported"
                    )
                container.append(parse_scalar(item))
                position += 1
                continue

            if content.startswith("- ") or ":" not in content:
                raise ValidationInputError(
                    f"{path}:{line_number}: expected a mapping entry"
                )
            key, raw_value = content.split(":", 1)
            key = key.strip()
            if not re.fullmatch(r"[A-Za-z0-9_.-]+", key):
                raise ValidationInputError(
                    f"{path}:{line_number}: unsupported key syntax: {key}"
                )
            if key in container:
                raise ValidationInputError(
                    f"{path}:{line_number}: duplicate key at this level: {key}"
                )
            raw_value = raw_value.strip()
            position += 1
            if raw_value:
                container[key] = parse_scalar(raw_value)
            elif position < len(logical_lines) and logical_lines[position][1] > indent:
                child_indent = logical_lines[position][1]
                if child_indent != indent + 2:
                    child_line = logical_lines[position][0]
                    raise ValidationInputError(
                        f"{path}:{child_line}: nested blocks must indent by two spaces"
                    )
                container[key], position = parse_block(position, child_indent)
            else:
                container[key] = {}

        return container, position

    root_indent = logical_lines[0][1]
    if root_indent != 0:
        raise ValidationInputError(f"{path}: root mapping must start at column 1")
    parsed, final_position = parse_block(0, 0)
    if final_position != len(logical_lines):
        line_number = logical_lines[final_position][0]
        raise ValidationInputError(f"{path}:{line_number}: unparsed content remains")
    if not isinstance(parsed, dict):
        raise ValidationInputError(f"{path}: root value must be a mapping")
    return parsed


def get_nested(config: dict[str, Any], path: str) -> Any:
    current: Any = config
    for key in path.split("."):
        if not isinstance(current, dict) or key not in current:
            raise KeyError(path)
        current = current[key]
    return current


def expect_value(
    config: dict[str, Any], path: str, expected: Any, reporter: Reporter
) -> None:
    try:
        observed = get_nested(config, path)
    except KeyError:
        reporter.fail("CONFIG_MISSING_KEY", f"Required setting is absent: {path}")
        return
    if observed != expected:
        reporter.fail(
            "CONFIG_VALUE",
            f"{path} is {observed!r}; expected current-state value {expected!r}",
        )


def expect_number_range(
    config: dict[str, Any], path: str, minimum: float, maximum: float,
    reporter: Reporter,
) -> None:
    try:
        value = get_nested(config, path)
    except KeyError:
        reporter.fail("CONFIG_MISSING_KEY", f"Required setting is absent: {path}")
        return
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        reporter.fail("CONFIG_TYPE", f"{path} must be numeric, not {value!r}")
    elif not minimum <= value <= maximum:
        reporter.fail(
            "CONFIG_RANGE", f"{path}={value} is outside [{minimum}, {maximum}]"
        )


def validate_config_schema(config: dict[str, Any], reporter: Reporter) -> None:
    expect_value(config, "schema_version", 1, reporter)
    expect_value(config, "snapshot.status", "documentation_snapshot", reporter)
    expect_value(config, "snapshot.scripts_read_this_file", False, reporter)
    expect_value(config, "study.sample_ids", EXPECTED_SAMPLES, reporter)
    expect_value(config, "study.structural_domains", EXPECTED_DOMAINS, reporter)
    expect_value(
        config,
        "public_data_archive.record_url",
        "https://zenodo.org/records/22058284",
        reporter,
    )
    expect_value(
        config, "public_data_archive.doi", "10.5281/zenodo.22058284", reporter
    )
    expect_value(
        config,
        "public_data_archive.expected_download_directory",
        "data/processed",
        reporter,
    )
    expect_value(
        config, "public_data_archive.file_groups.space_ranger_web_summary_html", 14, reporter
    )
    expect_value(
        config, "public_data_archive.file_groups.loupe_browser_cloupe", 14, reporter
    )
    expect_value(
        config, "public_data_archive.file_groups.velocyto_loom", 14, reporter
    )
    expect_value(
        config, "public_data_archive.file_groups.individual_visium_seurat_rds", 14, reporter
    )
    expect_value(
        config, "public_data_archive.complete_space_ranger_outs_included", False, reporter
    )
    expect_value(
        config, "scina.source_sheet", "Supplementary Table 9", reporter
    )
    expect_value(config, "scina.source_column_mapping.gene_id", "gene", reporter)
    expect_value(
        config, "scina.source_column_mapping.cell_type", "clusterName", reporter
    )
    expect_value(config, "scina.source_has_explicit_marker_rank", False, reporter)

    try:
        suffixes = get_nested(config, "study.barcode_suffixes")
    except KeyError:
        reporter.fail("CONFIG_MISSING_KEY", "study.barcode_suffixes is absent")
        suffixes = {}
    expected_suffixes = {sample: index for index, sample in enumerate(EXPECTED_SAMPLES, 1)}
    if suffixes != expected_suffixes:
        reporter.fail(
            "BARCODE_SUFFIX_CONFIG",
            "Configured sample-to-barcode suffix mapping is not the fixed 1-14 mapping",
        )

    ranges = [
        ("visium_qc.gene_filter.minimum_total_reads_across_retained_spots", 1, 10**9),
        ("visium_normalization_and_integration.pca.dimensions_tested", 2, 10000),
        ("visium_normalization_and_integration.pca.dimensions_used", 2, 10000),
        ("marker_analysis.min_detection_fraction", 0, 1),
        ("marker_analysis.log_fold_change_threshold", 0, 100),
        ("marker_analysis.adjusted_p_value_cutoff", 0, 1),
        ("developmental_trends_and_go.trend_clusters", 2, 100),
        ("developmental_trends_and_go.go_genes_per_cluster", 1, 100000),
        ("rna_velocity.min_shared_counts", 1, 10**9),
        ("rna_velocity.top_genes", 1, 10**7),
        ("rna_velocity.pca_dimensions", 2, 10000),
        ("rna_velocity.neighbors", 2, 10000),
        ("monocle3.minimum_total_umi_per_spot", 1, 10**12),
        ("monocle3.umap_neighbors", 2, 10000),
        ("monocle3.umap_min_dist", 0, 1),
        ("scrna_reference.cluster_resolution", 0, 100),
        ("spotlight_mapping.minimum_reported_proportion", 0, 1),
        ("spotlight_mapping.high_purity_cutoff", 0, 1),
    ]
    for path, minimum, maximum in ranges:
        expect_number_range(config, path, minimum, maximum, reporter)

    try:
        tested = get_nested(
            config, "visium_normalization_and_integration.pca.dimensions_tested"
        )
        used = get_nested(
            config, "visium_normalization_and_integration.pca.dimensions_used"
        )
        if isinstance(tested, int) and isinstance(used, int) and used > tested:
            reporter.fail(
                "PCA_DIMENSIONS", "PCA dimensions used cannot exceed dimensions tested"
            )
    except KeyError:
        pass

    try:
        config_paths = get_nested(config, "paths")
    except KeyError:
        reporter.fail("CONFIG_MISSING_KEY", "paths mapping is absent")
        config_paths = {}
    for name, value in config_paths.items():
        if not isinstance(value, str):
            reporter.fail("PATH_TYPE", f"paths.{name} must be a string")
            continue
        if re.match(r"^[A-Za-z]:[/\\]", value) or value.startswith(("/", "\\")):
            reporter.fail("ABSOLUTE_PATH", f"paths.{name} contains an absolute path")
        normalized = value.replace("\\", "/")
        if ".." in PurePosixPath(normalized).parts:
            reporter.fail("PARENT_PATH", f"paths.{name} escapes the project root")

    if reporter.count("FAIL") == 0:
        reporter.passed(
            "CONFIG_SCHEMA",
            "Configuration keys, fixed sample mapping, value ranges, and portable paths are valid",
        )


def validate_sample_manifest(
    project_root: Path, config: dict[str, Any], reporter: Reporter
) -> None:
    path = project_root / "config" / "spaceranger_samples.csv"
    if not path.is_file():
        reporter.fail("SAMPLE_MANIFEST", f"Sample manifest is missing: {path}")
        return
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {
        "sample_id", "barcode_suffix", "spaceranger_sample_name", "slide_id",
        "capture_area_from_script", "image_region_in_filename", "local_cores",
        "local_memory_gb", "input_paths_status", "verification_note",
    }
    missing = required.difference(rows[0].keys() if rows else set())
    if missing:
        reporter.fail(
            "SAMPLE_MANIFEST_COLUMNS",
            "Sample manifest is missing column(s): " + ", ".join(sorted(missing)),
        )
        return
    observed_samples = [row["sample_id"] for row in rows]
    configured_samples = get_nested(config, "study.sample_ids")
    if observed_samples != configured_samples:
        reporter.fail(
            "SAMPLE_ORDER", "Sample manifest order does not match study.sample_ids"
        )
    if len(set(observed_samples)) != len(observed_samples):
        reporter.fail("SAMPLE_DUPLICATE", "Sample manifest contains duplicate sample IDs")

    bad_suffixes: list[str] = []
    area_mismatches: list[str] = []
    for index, row in enumerate(rows, 1):
        try:
            suffix = int(row["barcode_suffix"])
        except ValueError:
            bad_suffixes.append(row["sample_id"])
            continue
        if suffix != index:
            bad_suffixes.append(row["sample_id"])
        if row["capture_area_from_script"] != row["image_region_in_filename"]:
            area_mismatches.append(row["sample_id"])
    if bad_suffixes:
        reporter.fail(
            "SAMPLE_SUFFIX", "Invalid barcode suffixes for: " + ", ".join(bad_suffixes)
        )
    if area_mismatches:
        reporter.warn(
            "CAPTURE_AREA_MISMATCH",
            "Legacy command and image filename disagree for: "
            + ", ".join(area_mismatches),
        )
    if rows and not missing and not bad_suffixes and observed_samples == configured_samples:
        reporter.passed(
            "SAMPLE_MANIFEST", "14-sample order, uniqueness, and barcode suffixes are valid"
        )


def validate_metadata(
    project_root: Path, config: dict[str, Any], reporter: Reporter
) -> None:
    path = project_root / "data" / "metadata" / "metadata.csv"
    if not path.is_file():
        reporter.fail("METADATA_FILE", f"Metadata file is missing: {path}")
        return

    samples = get_nested(config, "study.sample_ids")
    domains = get_nested(config, "study.structural_domains")
    suffixes = get_nested(config, "study.barcode_suffixes")
    observed_samples: Counter[str] = Counter()
    observed_domains: Counter[str] = Counter()
    seen_barcodes: set[str] = set()
    duplicate_examples: list[str] = []
    invalid_samples: set[str] = set()
    invalid_domains: set[str] = set()
    suffix_examples: list[str] = []
    section_examples: list[str] = []
    blank_examples: list[str] = []
    row_count = 0

    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing_columns = REQUIRED_METADATA_COLUMNS.difference(fields)
        if missing_columns:
            reporter.fail(
                "METADATA_COLUMNS",
                "metadata.csv is missing column(s): "
                + ", ".join(sorted(missing_columns)),
            )
            return
        for row_number, row in enumerate(reader, start=2):
            row_count += 1
            barcode = (row.get("Barcode") or "").strip()
            sample = (row.get("sample_id") or "").strip()
            domain = (row.get("domains") or "").strip()
            section_id = (row.get("section_id") or "").strip()
            cluster = (row.get("harmony_clusters") or "").strip()
            if not all((barcode, sample, domain, section_id, cluster)):
                if len(blank_examples) < 5:
                    blank_examples.append(str(row_number))
                continue
            if barcode in seen_barcodes and len(duplicate_examples) < 5:
                duplicate_examples.append(barcode)
            seen_barcodes.add(barcode)
            observed_samples[sample] += 1
            observed_domains[domain] += 1
            if sample not in samples:
                invalid_samples.add(sample)
            if domain not in domains:
                invalid_domains.add(domain)
            if sample in suffixes and not barcode.endswith(f"_1_{suffixes[sample]}"):
                if len(suffix_examples) < 5:
                    suffix_examples.append(barcode)
            if not section_id.startswith(f"{sample}_") and len(section_examples) < 5:
                section_examples.append(section_id)

    if row_count == 0:
        reporter.fail("METADATA_EMPTY", "metadata.csv contains no data rows")
    if blank_examples:
        reporter.fail(
            "METADATA_BLANK", "Required metadata values are blank at row(s): "
            + ", ".join(blank_examples)
        )
    if duplicate_examples:
        reporter.fail(
            "METADATA_DUPLICATE_BARCODE",
            "Duplicate barcode example(s): " + ", ".join(duplicate_examples),
        )
    if invalid_samples:
        reporter.fail(
            "METADATA_SAMPLE", "Unexpected sample ID(s): " + ", ".join(sorted(invalid_samples))
        )
    if invalid_domains:
        reporter.fail(
            "METADATA_DOMAIN", "Unexpected domain(s): " + ", ".join(sorted(invalid_domains))
        )
    missing_samples = [sample for sample in samples if observed_samples[sample] == 0]
    missing_domains = [domain for domain in domains if observed_domains[domain] == 0]
    if missing_samples:
        reporter.fail("METADATA_SAMPLE_COVERAGE", "No rows for: " + ", ".join(missing_samples))
    if missing_domains:
        reporter.fail("METADATA_DOMAIN_COVERAGE", "No rows for: " + ", ".join(missing_domains))
    if suffix_examples:
        reporter.fail(
            "METADATA_BARCODE_SUFFIX",
            "Barcode suffix does not match sample mapping; example(s): "
            + ", ".join(suffix_examples),
        )
    if section_examples:
        reporter.fail(
            "METADATA_SECTION_ID",
            "section_id does not begin with sample_id; example(s): "
            + ", ".join(section_examples),
        )

    metadata_failure_codes = {
        "METADATA_EMPTY", "METADATA_BLANK", "METADATA_DUPLICATE_BARCODE",
        "METADATA_SAMPLE", "METADATA_DOMAIN", "METADATA_SAMPLE_COVERAGE",
        "METADATA_DOMAIN_COVERAGE", "METADATA_BARCODE_SUFFIX", "METADATA_SECTION_ID",
    }
    if not any(
        item.severity == "FAIL" and item.code in metadata_failure_codes
        for item in reporter.findings
    ):
        reporter.passed(
            "METADATA_INTEGRITY",
            f"{row_count:,} rows have valid required fields, unique barcodes, samples, domains, and suffixes",
        )


def validate_present_files(
    project_root: Path, config: dict[str, Any], reporter: Reporter
) -> None:
    required_files = [
        "data/metadata/metadata.csv",
        "data/reference/maize_mitochondrial_genes.txt",
        "data/reference/maize_plastid_genes.txt",
        "data/processed/maize_shoot_14samples_SCT_harmony_seurat_v5.rds",
        "data/processed/sc_merged_filter_SCT2_inte.rds",
        "data/reference/developmental_trends/sub_gene.csv",
        "data/reference/developmental_trends/zea_go2.csv",
        "data/reference/developmental_trends/maize_id_name.csv",
    ]
    samples = get_nested(config, "study.sample_ids")
    required_files.extend(
        f"data/processed/{sample}_seurat_v5.rds" for sample in samples
    )
    missing: list[str] = []
    empty: list[str] = []
    for relative in required_files:
        path = project_root / relative
        if not path.is_file():
            missing.append(relative)
        elif path.stat().st_size == 0:
            empty.append(relative)
    if missing:
        reporter.fail("PRESENT_FILE_MISSING", "Declared-present file(s) missing: " + ", ".join(missing))
    if empty:
        reporter.fail("PRESENT_FILE_EMPTY", "Required file(s) are empty: " + ", ".join(empty))
    if not missing and not empty:
        reporter.passed(
            "PRESENT_FILES",
            f"All {len(required_files)} files declared present exist and are non-empty",
        )


def regex_assignment(name: str, value: Any) -> str:
    if isinstance(value, bool):
        rendered = "TRUE" if value else "FALSE"
    elif isinstance(value, str):
        rendered = rf'["\']{re.escape(value)}["\']'
    else:
        rendered = re.escape(str(value)) + "L?"
    return rf"\b{re.escape(name)}\s*<-\s*{rendered}(?=\s|$)"


def validate_script_consistency(
    project_root: Path, config: dict[str, Any], reporter: Reporter
) -> None:
    checks: list[tuple[str, str, str]] = []

    integration = "scripts/R/04_merge_integration/01_merge_14_samples_SCT_Harmony_Seurat_v5.R"
    checks.extend([
        (integration, "max_pcs_to_test", regex_assignment("max_pcs_to_test", get_nested(config, "visium_normalization_and_integration.pca.dimensions_tested"))),
        (integration, "n_pcs_use", regex_assignment("n_pcs_use", get_nested(config, "visium_normalization_and_integration.pca.dimensions_used"))),
        (integration, "minimum_total_reads_per_gene", regex_assignment("minimum_total_reads_per_gene", get_nested(config, "visium_qc.gene_filter.minimum_total_reads_across_retained_spots"))),
        (integration, "random_seed", regex_assignment("random_seed", get_nested(config, "visium_normalization_and_integration.random_seed"))),
        (integration, "UMAP n.neighbors", r"n\.neighbors\s*=\s*30L"),
        (integration, "UMAP min.dist", r"min\.dist\s*=\s*0\.3"),
        (integration, "UMAP metric", r'metric\s*=\s*["\']cosine["\']'),
    ])

    marker = "scripts/R/05_marker_analysis/01_find_markers_harmony_clusters_SCT_Seurat_v5.R"
    checks.extend([
        (marker, "marker_min_pct", regex_assignment("marker_min_pct", get_nested(config, "marker_analysis.min_detection_fraction"))),
        (marker, "marker_logfc_threshold", regex_assignment("marker_logfc_threshold", get_nested(config, "marker_analysis.log_fold_change_threshold"))),
        (marker, "marker_adjusted_p_cutoff", regex_assignment("marker_adjusted_p_cutoff", get_nested(config, "marker_analysis.adjusted_p_value_cutoff"))),
    ])

    trends = "scripts/R/07_developmental_trends_GO/01_developmental_expression_trends_and_GO_Figure_9_Seurat_v5.R"
    checks.extend([
        (trends, "number_of_clusters", regex_assignment("number_of_clusters", get_nested(config, "developmental_trends_and_go.trend_clusters"))),
        (trends, "representatives_per_cluster", regex_assignment("representatives_per_cluster", get_nested(config, "developmental_trends_and_go.representative_genes_per_cluster"))),
        (trends, "go_genes_per_cluster", regex_assignment("go_genes_per_cluster", get_nested(config, "developmental_trends_and_go.go_genes_per_cluster"))),
        (
            trends,
            "go_pvalue_cutoff",
            r"\bgo_pvalue_cutoff\s*<-\s*1[eE]-?0*6(?=\s|$)",
        ),
    ])

    monocle = "scripts/R/09_monocle3_pseudotime/01_monocle3_pseudotime_14_samples_Figure_10C.R"
    checks.extend([
        (monocle, "umi_cutoff", regex_assignment("umi_cutoff", get_nested(config, "monocle3.minimum_total_umi_per_spot"))),
        (monocle, "num_dim", regex_assignment("num_dim", get_nested(config, "monocle3.pca_dimensions"))),
        (monocle, "alignment_column", regex_assignment("alignment_column", get_nested(config, "monocle3.alignment_group"))),
        (monocle, "umap_n_neighbors", regex_assignment("umap_n_neighbors", get_nested(config, "monocle3.umap_neighbors"))),
        (monocle, "umap_min_dist", regex_assignment("umap_min_dist", get_nested(config, "monocle3.umap_min_dist"))),
        (monocle, "random_seed", regex_assignment("random_seed", get_nested(config, "monocle3.random_seed"))),
    ])

    scrna = "scripts/R/10_scRNA_reference_integration/01_prepare_maize_scRNA_reference_SCT_Harmony_Seurat_v5.R"
    checks.extend([
        (scrna, "min_cells", regex_assignment("min_cells", get_nested(config, "scrna_reference.create_seurat_object.min_cells"))),
        (scrna, "min_features", regex_assignment("min_features", get_nested(config, "scrna_reference.create_seurat_object.min_features"))),
        (scrna, "variable_features", regex_assignment("variable_features", get_nested(config, "scrna_reference.sctransform.variable_features"))),
        (scrna, "integration_pcs", regex_assignment("integration_pcs", get_nested(config, "scrna_reference.harmony_dimensions"))),
        (scrna, "cluster_resolution", regex_assignment("cluster_resolution", get_nested(config, "scrna_reference.cluster_resolution"))),
        (scrna, "random_seed", regex_assignment("random_seed", get_nested(config, "scrna_reference.random_seed"))),
    ])

    failures: list[str] = []
    cache: dict[str, str] = {}
    for relative, label, pattern in checks:
        path = project_root / relative
        if not path.is_file():
            failures.append(f"{relative}: missing")
            continue
        text = cache.setdefault(relative, path.read_text(encoding="utf-8-sig"))
        if not re.search(pattern, text):
            failures.append(f"{relative}: {label}")
    if failures:
        reporter.fail(
            "SCRIPT_CONFIG_DRIFT",
            "Configured value not found in current script assignment(s): "
            + "; ".join(failures),
        )
    else:
        reporter.passed(
            "SCRIPT_CONFIG_SNAPSHOT",
            f"{len(checks)} important R-script assignments match the configuration snapshot",
        )

    integration_text = cache.get(integration, "")
    repeated_umap_settings = {
        "n.neighbors = 30L": r"n\.neighbors\s*=\s*30L",
        "min.dist = 0.3": r"min\.dist\s*=\s*0\.3",
        'metric = "cosine"': r'metric\s*=\s*["\']cosine["\']',
    }
    incomplete_umap = [
        label for label, pattern in repeated_umap_settings.items()
        if len(re.findall(pattern, integration_text)) != 2
    ]
    if incomplete_umap:
        reporter.fail(
            "UMAP_EXPLICIT_SETTINGS",
            "Each pre/post-integration RunUMAP call must explicitly set: "
            + ", ".join(incomplete_umap),
        )
    else:
        reporter.passed(
            "UMAP_EXPLICIT_SETTINGS",
            "Both pre- and post-integration UMAP calls explicitly record neighbors, minimum distance, and metric",
        )

    velocity_path = project_root / "scripts/python/08_RNA_velocity/01_scvelo_dynamical_RNA_velocity.py"
    if not velocity_path.is_file():
        reporter.fail("VELOCITY_SCRIPT", "RNA-velocity script is missing")
    else:
        velocity_text = velocity_path.read_text(encoding="utf-8-sig")
        defaults = {
            "min-shared-counts": get_nested(config, "rna_velocity.min_shared_counts"),
            "n-top-genes": get_nested(config, "rna_velocity.top_genes"),
            "n-pcs": get_nested(config, "rna_velocity.pca_dimensions"),
            "n-neighbors": get_nested(config, "rna_velocity.neighbors"),
            "max-iter": get_nested(config, "rna_velocity.recover_dynamics_max_iterations"),
            "n-jobs": get_nested(config, "rna_velocity.parallel_jobs"),
            "diagnostic-unspliced-threshold": get_nested(
                config, "rna_velocity.low_unspliced_fraction_diagnostic_threshold"
            ),
        }
        missing_defaults = []
        for option, value in defaults.items():
            pattern = (
                rf'add_argument\(\s*["\']--{re.escape(option)}["\'].*?'
                rf'default\s*=\s*{re.escape(str(value))}\b'
            )
            if not re.search(pattern, velocity_text, flags=re.DOTALL):
                missing_defaults.append(option)
        if missing_defaults:
            reporter.fail(
                "VELOCITY_CONFIG_DRIFT",
                "RNA-velocity CLI defaults differ or could not be located: "
                + ", ".join(missing_defaults),
            )
        else:
            reporter.passed(
                "VELOCITY_CONFIG_SNAPSHOT",
                "RNA-velocity command-line defaults match the configuration snapshot",
            )


def report_unresolved_settings(config: dict[str, Any], reporter: Reporter) -> None:
    unresolved_paths = [
        "visium_qc.spot_level_thresholds.detected_genes_min",
        "visium_qc.spot_level_thresholds.detected_genes_max",
        "visium_qc.spot_level_thresholds.total_umi_min",
        "visium_qc.spot_level_thresholds.total_umi_max",
        "visium_qc.spot_level_thresholds.mitochondrial_percent_max",
        "visium_qc.spot_level_thresholds.chloroplast_percent_max",
        "visium_normalization_and_integration.harmony.grouping_variable",
        "visium_normalization_and_integration.clustering.resolution",
        "space_ranger_reference.fasta",
        "space_ranger_reference.input_gtf",
    ]
    unresolved = []
    for path in unresolved_paths:
        try:
            if get_nested(config, path) is None:
                unresolved.append(path)
        except KeyError:
            unresolved.append(path)
    if unresolved:
        reporter.warn(
            "UNRESOLVED_SETTINGS",
            f"{len(unresolved)} publication-relevant setting(s) remain null: "
            + ", ".join(unresolved),
        )

    try:
        missing_blocking = get_nested(config, "data_readiness.missing_blocking")
    except KeyError:
        missing_blocking = []
    if missing_blocking:
        reporter.warn(
            "MISSING_BLOCKING_INPUTS",
            "End-to-end rebuild inputs remain unavailable: "
            + ", ".join(str(item) for item in missing_blocking),
        )
    reporter.warn(
        "CONFIG_NOT_AUTHORITATIVE",
        "Scripts still contain their own parameter assignments and do not load current_state.yaml",
    )


def write_json_report(path: Path, reporter: Reporter) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "summary": {
            "passed": reporter.count("PASS"),
            "warnings": reporter.count("WARN"),
            "failures": reporter.count("FAIL"),
        },
        "findings": [asdict(item) for item in reporter.findings],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def run_validation(project_root: Path) -> tuple[dict[str, Any], Reporter]:
    reporter = Reporter()
    config_path = project_root / "config" / "current_state.yaml"
    try:
        config = load_yaml_subset(config_path)
    except (OSError, ValidationInputError) as error:
        reporter.fail("CONFIG_PARSE", str(error))
        return {}, reporter
    validate_config_schema(config, reporter)
    validate_sample_manifest(project_root, config, reporter)
    validate_metadata(project_root, config, reporter)
    validate_present_files(project_root, config, reporter)
    validate_script_consistency(project_root, config, reporter)
    report_unresolved_settings(config, reporter)
    return config, reporter


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Validate configuration, metadata, declared inputs, and script consistency."
    )
    parser.add_argument(
        "--project-root", type=Path, default=project_root_from_script(),
        help="Repository root; defaults to the root containing this script.",
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Treat warnings (including unresolved scientific settings) as a nonzero result.",
    )
    parser.add_argument(
        "--json-report", type=Path, default=None,
        help="Optional path for a machine-readable validation report.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_arguments()
    project_root = args.project_root.resolve()
    _, reporter = run_validation(project_root)
    reporter.print_report()
    if args.json_report is not None:
        report_path = args.json_report
        if not report_path.is_absolute():
            report_path = project_root / report_path
        write_json_report(report_path, reporter)
        print(f"JSON report: {report_path}")
    if reporter.count("FAIL"):
        return 1
    if args.strict and reporter.count("WARN"):
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
