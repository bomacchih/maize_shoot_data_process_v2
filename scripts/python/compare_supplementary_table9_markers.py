#!/usr/bin/env python3
"""Compare the curated SCINA marker CSV with Supplementary Table 9.

The XLSX reader uses only the Python standard library so this provenance check
does not add a runtime dependency to the biological analysis environment.
"""

from __future__ import annotations

import argparse
import csv
import json
import posixpath
import re
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple
from xml.etree import ElementTree


SHEET_NAME = "Supplementary Table 9"
WORKBOOK_XML = "xl/workbook.xml"
WORKBOOK_RELS_XML = "xl/_rels/workbook.xml.rels"
RELATIONSHIP_ID = (
    "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
)


def local_name(tag: str) -> str:
    return tag.rsplit("}", 1)[-1]


def normalize_cell_type(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9]+", "_", value.strip()).strip("_")


def column_index(cell_reference: str) -> int:
    letters = re.match(r"[A-Za-z]+", cell_reference)
    if letters is None:
        return -1
    value = 0
    for letter in letters.group(0).upper():
        value = value * 26 + ord(letter) - ord("A") + 1
    return value - 1


def shared_strings(archive: zipfile.ZipFile) -> List[str]:
    path = "xl/sharedStrings.xml"
    if path not in archive.namelist():
        return []
    values: List[str] = []
    with archive.open(path) as handle:
        for _, element in ElementTree.iterparse(handle, events=("end",)):
            if local_name(element.tag) != "si":
                continue
            values.append(
                "".join(
                    node.text or ""
                    for node in element.iter()
                    if local_name(node.tag) == "t"
                )
            )
            element.clear()
    return values


def sheet_xml_path(archive: zipfile.ZipFile, sheet_name: str) -> str:
    workbook = ElementTree.parse(archive.open(WORKBOOK_XML)).getroot()
    relationship_id = None
    for element in workbook.iter():
        if local_name(element.tag) == "sheet" and element.get("name") == sheet_name:
            relationship_id = element.get(RELATIONSHIP_ID)
            break
    if relationship_id is None:
        raise ValueError(f"Workbook has no sheet named {sheet_name!r}.")

    relationships = ElementTree.parse(
        archive.open(WORKBOOK_RELS_XML)
    ).getroot()
    target = None
    for element in relationships:
        if element.get("Id") == relationship_id:
            target = element.get("Target")
            break
    if target is None:
        raise ValueError(f"Could not resolve worksheet relationship {relationship_id!r}.")

    if target.startswith("/"):
        return target.lstrip("/")
    return posixpath.normpath(posixpath.join("xl", target))


def cell_text(cell: ElementTree.Element, strings: Sequence[str]) -> str:
    cell_type = cell.get("t")
    if cell_type == "inlineStr":
        return "".join(
            node.text or ""
            for node in cell.iter()
            if local_name(node.tag) == "t"
        ).strip()

    value_node = next(
        (node for node in cell if local_name(node.tag) == "v"), None
    )
    if value_node is None or value_node.text is None:
        return ""
    raw_value = value_node.text.strip()
    if cell_type == "s":
        return strings[int(raw_value)].strip()
    return raw_value


def read_table9_rows(workbook_path: Path) -> List[Dict[str, object]]:
    with zipfile.ZipFile(workbook_path) as archive:
        strings = shared_strings(archive)
        worksheet_path = sheet_xml_path(archive, SHEET_NAME)
        rows: List[Dict[str, object]] = []
        last_nonblank_row = 0
        with archive.open(worksheet_path) as handle:
            for _, element in ElementTree.iterparse(handle, events=("end",)):
                if local_name(element.tag) != "row":
                    continue
                row_number = int(element.get("r", "0"))
                if row_number < 4:
                    element.clear()
                    continue
                if last_nonblank_row and row_number - last_nonblank_row >= 1000:
                    break

                values = ["", "", ""]
                for cell in element:
                    if local_name(cell.tag) != "c":
                        continue
                    index = column_index(cell.get("r", ""))
                    if 0 <= index < 3:
                        values[index] = cell_text(cell, strings)

                if any(values):
                    last_nonblank_row = row_number
                    gene_id, source_cell_type, celltype_id = values
                    rows.append(
                        {
                            "source_row": row_number,
                            "gene_id": gene_id,
                            "cell_type": normalize_cell_type(source_cell_type),
                            "source_cell_type": source_cell_type,
                            "celltype_id": celltype_id,
                        }
                    )
                element.clear()
    return rows


def read_marker_rows(marker_path: Path) -> List[Dict[str, str]]:
    with marker_path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    required = {"gene_id", "cell_type", "marker_rank", "source_cell_type"}
    missing = required.difference(rows[0] if rows else {})
    if missing:
        raise ValueError(
            "Marker CSV is missing required column(s): " + ", ".join(sorted(missing))
        )
    return rows


def unique_in_order(rows: Iterable[Dict[str, object]]) -> List[Dict[str, object]]:
    seen = set()
    result = []
    for row in rows:
        pair = (row["gene_id"], row["cell_type"])
        if pair not in seen:
            seen.add(pair)
            result.append(row)
    return result


def compare(table_rows: List[Dict[str, object]], marker_rows: List[Dict[str, str]]):
    malformed = [
        row for row in table_rows if not row["gene_id"] or not row["cell_type"]
    ]
    valid_rows = [row for row in table_rows if row not in malformed]
    unique_pairs = unique_in_order(valid_rows)

    gene_to_types: Dict[object, set] = defaultdict(set)
    for row in unique_pairs:
        gene_to_types[row["gene_id"]].add(row["cell_type"])
    shared_genes = {
        gene_id for gene_id, cell_types in gene_to_types.items()
        if len(cell_types) > 1
    }
    exclusive_pairs = [
        row for row in unique_pairs if row["gene_id"] not in shared_genes
    ]

    marker_pairs = {(row["gene_id"], row["cell_type"]) for row in marker_rows}
    unique_pair_set = {
        (row["gene_id"], row["cell_type"]) for row in unique_pairs
    }
    exclusive_pair_set = {
        (row["gene_id"], row["cell_type"]) for row in exclusive_pairs
    }

    marker_gene_to_types: Dict[str, set] = defaultdict(set)
    for row in marker_rows:
        marker_gene_to_types[row["gene_id"]].add(row["cell_type"])

    source_sequences: Dict[object, list] = defaultdict(list)
    for row in exclusive_pairs:
        pair = (row["gene_id"], row["cell_type"])
        if pair in marker_pairs:
            source_sequences[row["cell_type"]].append(row["gene_id"])

    marker_sequences: Dict[str, list] = defaultdict(list)
    for row in sorted(
        marker_rows, key=lambda item: (item["cell_type"], int(item["marker_rank"]))
    ):
        marker_sequences[row["cell_type"]].append(row["gene_id"])

    order_matches = {
        cell_type: marker_sequences[cell_type] == source_sequences[cell_type]
        for cell_type in sorted(set(marker_sequences) | set(source_sequences))
    }
    comparison_passed = all(
        (
            not malformed,
            len(marker_pairs) == len(marker_rows),
            not any(len(types) > 1 for types in marker_gene_to_types.values()),
            marker_pairs <= exclusive_pair_set,
            all(order_matches.values()),
        )
    )

    return {
        "comparison_passed": comparison_passed,
        "source_sheet": SHEET_NAME,
        "table9_raw_rows": len(table_rows),
        "table9_malformed_rows": len(malformed),
        "table9_duplicate_pair_rows": len(valid_rows) - len(unique_pairs),
        "table9_unique_gene_celltype_pairs": len(unique_pairs),
        "table9_unique_genes": len(gene_to_types),
        "table9_genes_shared_across_cell_types": len(shared_genes),
        "table9_exclusive_genes": len(exclusive_pair_set),
        "marker_rows": len(marker_rows),
        "marker_cell_types": len(marker_sequences),
        "marker_pairs_absent_from_table9": len(marker_pairs - unique_pair_set),
        "marker_pairs_not_celltype_exclusive_in_table9": len(
            marker_pairs - exclusive_pair_set
        ),
        "table9_exclusive_pairs_not_in_marker_file": len(
            exclusive_pair_set - marker_pairs
        ),
        "marker_source_order_preserved": all(order_matches.values()),
        "marker_counts_by_cell_type": dict(
            sorted(Counter(row["cell_type"] for row in marker_rows).items())
        ),
    }


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path)
    parser.add_argument("marker_csv", type=Path)
    parser.add_argument("--json-report", type=Path)
    return parser.parse_args()


def main() -> int:
    arguments = parse_arguments()
    table_rows = read_table9_rows(arguments.workbook)
    marker_rows = read_marker_rows(arguments.marker_csv)
    result = compare(table_rows, marker_rows)
    output = json.dumps(result, indent=2, sort_keys=True)
    print(output)
    if arguments.json_report:
        arguments.json_report.parent.mkdir(parents=True, exist_ok=True)
        arguments.json_report.write_text(output + "\n", encoding="utf-8")
    return 0 if result["comparison_passed"] else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, zipfile.BadZipFile) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(2)
