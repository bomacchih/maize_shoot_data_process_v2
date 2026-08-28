#!/usr/bin/env python3
"""Audit a Git repository before changing its GitHub visibility to public.

The audit scans the current tracked tree and reachable Git history. Potential
secret values are never printed or written to reports; only finding types and
file/line locations are recorded.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from collections import defaultdict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


MIB = 1024 * 1024
AUDIT_SCRIPT = "scripts/python/validate_public_readiness.py"
SCAN_CHUNK_CHARACTERS = 16_384
MAX_LOCATIONS_PER_CODE = 200

TEXT_EXTENSIONS = {
    ".cff", ".cfg", ".css", ".csv", ".html", ".ini", ".js", ".json",
    ".md", ".py", ".r", ".rmd", ".rst", ".sh", ".toml", ".tsv",
    ".txt", ".xml", ".yaml", ".yml",
}
TEXT_BASENAMES = {
    ".gitattributes", ".gitignore", "citation", "license", "licence",
    "makefile", "readme", "security",
}

SENSITIVE_DATA_SUFFIXES = (
    ".rds", ".rdata", ".h5", ".h5ad", ".loom", ".bam", ".bai",
    ".fastq", ".fastq.gz", ".fq", ".fq.gz", ".cloupe",
)
SENSITIVE_CREDENTIAL_SUFFIXES = (
    ".env", ".pem", ".p12", ".pfx", ".key", ".kdbx",
)

SECRET_PATTERNS = {
    "PRIVATE_KEY": re.compile(
        r"-----BEGIN (?:RSA |EC |DSA |OPENSSH )?PRIVATE KEY-----"
    ),
    "GITHUB_TOKEN": re.compile(
        r"\b(?:gh[pousr]_[A-Za-z0-9]{30,}|github_pat_[A-Za-z0-9_]{50,})\b"
    ),
    "AWS_ACCESS_KEY": re.compile(r"\b(?:AKIA|ASIA)[A-Z0-9]{16}\b"),
    "SLACK_TOKEN": re.compile(r"\bxox[baprs]-[A-Za-z0-9-]{20,}\b"),
    "OPENAI_KEY": re.compile(r"\bsk-(?:proj-)?[A-Za-z0-9_-]{20,}\b"),
    "GENERIC_SECRET_ASSIGNMENT": re.compile(
        r"(?i)\b(?:password|passwd|api[_-]?key|access[_-]?token|secret)\b"
        r"\s*[:=]\s*[\"']?[A-Za-z0-9+/=_-]{12,}"
    ),
}

WINDOWS_PERSONAL_PATH = re.compile(
    r"(?i)\b[A-Z]:[\\/](?:Users|Documents and Settings)[\\/][^\\/\s\"']+"
)
UNIX_PERSONAL_PATH = re.compile(r"/(?:home|Users)/([^/\s\"']+)")
EMAIL_ADDRESS = re.compile(
    r"\b[A-Za-z0-9.!#$%&'*+/=?^_`{|}~-]+@"
    r"[A-Za-z0-9](?:[A-Za-z0-9-]{0,61}[A-Za-z0-9])?"
    r"(?:\.[A-Za-z0-9](?:[A-Za-z0-9-]{0,61}[A-Za-z0-9])?)+\b"
)
R_SLOT_EMAIL_LOOKALIKES = ("@meta.data", "@data.frame")
URL = re.compile(r"https?://[^\s<>()\"']+")

PLACEHOLDER_TERMS = {
    "changeme", "example", "placeholder", "replace_me", "token_here",
    "your_api", "your_key", "your_token",
}
GENERIC_UNIX_USERS = {"jdoe", "user", "username"}


class AuditError(RuntimeError):
    """Raised when the audit cannot inspect the repository safely."""


@dataclass
class Finding:
    severity: str
    code: str
    message: str
    locations: list[str] = field(default_factory=list)


class Reporter:
    def __init__(self) -> None:
        self.findings: list[Finding] = []

    def add(
        self,
        severity: str,
        code: str,
        message: str,
        locations: Iterable[str] = (),
    ) -> None:
        unique_locations = sorted(set(locations))
        self.findings.append(
            Finding(severity, code, message, unique_locations)
        )

    def count(self, severity: str) -> int:
        return sum(f.severity == severity for f in self.findings)

    def print_console(self) -> None:
        for finding in self.findings:
            print(f"[{finding.severity}] {finding.code}: {finding.message}")
            shown = finding.locations[:12]
            for location in shown:
                print(f"  - {location}")
            if len(finding.locations) > len(shown):
                print(f"  - ... {len(finding.locations) - len(shown)} more location(s)")
        print(
            "\nSummary: "
            f"{self.count('PASS')} passed, "
            f"{self.count('WARN')} warning(s), "
            f"{self.count('BLOCKER')} blocker(s)."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path.cwd())
    parser.add_argument(
        "--git-executable", default=os.environ.get("GIT_EXE", "git")
    )
    parser.add_argument("--json-report", type=Path)
    parser.add_argument("--markdown-report", type=Path)
    parser.add_argument(
        "--history-max-blob-mb",
        type=float,
        default=2.0,
        help="Maximum historical text-blob size scanned for secrets.",
    )
    parser.add_argument(
        "--skip-history", action="store_true",
        help="Skip full-history secret and sensitive-file scanning."
    )
    parser.add_argument(
        "--strict", action="store_true",
        help="Return nonzero for warnings as well as blockers."
    )
    return parser.parse_args()


def run_git(
    root: Path,
    git_executable: str,
    arguments: list[str],
    *,
    input_bytes: bytes | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[bytes]:
    command = [
        git_executable,
        "-c",
        f"safe.directory={root.as_posix()}",
        *arguments,
    ]
    try:
        result = subprocess.run(
            command,
            cwd=root,
            input=input_bytes,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
    except FileNotFoundError as exc:
        raise AuditError(
            f"Git executable not found: {git_executable}. "
            "Use --git-executable or set GIT_EXE."
        ) from exc
    if check and result.returncode != 0:
        error = result.stderr.decode("utf-8", errors="replace").strip()
        raise AuditError(f"Git command failed: {' '.join(arguments)}: {error}")
    return result


def decode_paths(payload: bytes) -> list[str]:
    return [
        item.decode("utf-8", errors="surrogateescape")
        for item in payload.split(b"\0")
        if item
    ]


def is_text_path(path: str) -> bool:
    name = Path(path).name.lower()
    suffix = Path(path).suffix.lower()
    return suffix in TEXT_EXTENSIONS or name in TEXT_BASENAMES


def has_sensitive_suffix(path: str) -> bool:
    lower = path.lower()
    return lower.endswith(SENSITIVE_DATA_SUFFIXES + SENSITIVE_CREDENTIAL_SUFFIXES)


def is_placeholder_line(line: str) -> bool:
    lower = line.lower()
    return any(term in lower for term in PLACEHOLDER_TERMS) or "<token>" in lower


def classify_text_findings(
    text: str,
    location_prefix: str,
) -> dict[str, set[str]]:
    """Return redacted finding locations for one text document."""
    found: dict[str, set[str]] = defaultdict(set)
    for line_number, line in enumerate(text.splitlines(), start=1):
        location = f"{location_prefix}:{line_number}"
        # Rendered HTML can contain multi-megabyte minified lines. Bounded
        # chunks avoid pathological regular-expression runtimes while keeping
        # the file/line location useful and secret values redacted.
        segments = (
            line[start : start + SCAN_CHUNK_CHARACTERS]
            for start in range(0, max(1, len(line)), SCAN_CHUNK_CHARACTERS)
        )
        for segment in segments:
            for code, pattern in SECRET_PATTERNS.items():
                if (
                    len(found.get(code, set())) < MAX_LOCATIONS_PER_CODE
                    and pattern.search(segment)
                ):
                    if (
                        code == "GENERIC_SECRET_ASSIGNMENT"
                        and is_placeholder_line(segment)
                    ):
                        continue
                    found[code].add(location)

            if (
                len(found.get("PERSONAL_ABSOLUTE_PATH", set())) < MAX_LOCATIONS_PER_CODE
                and WINDOWS_PERSONAL_PATH.search(segment)
            ):
                found["PERSONAL_ABSOLUTE_PATH"].add(location)
            for match in UNIX_PERSONAL_PATH.finditer(segment):
                if (
                    len(found.get("PERSONAL_ABSOLUTE_PATH", set())) < MAX_LOCATIONS_PER_CODE
                    and match.group(1).lower() not in GENERIC_UNIX_USERS
                ):
                    found["PERSONAL_ABSOLUTE_PATH"].add(location)

            email_matches = list(EMAIL_ADDRESS.finditer(segment))
            path_without_blob = location_prefix.split("@blob-", 1)[0].lower()
            is_r_source = path_without_blob.endswith((".r", ".rmd"))
            public_email_matches = [
                match for match in email_matches
                if not (
                    is_r_source
                    and any(
                        marker in match.group(0).lower()
                        for marker in R_SLOT_EMAIL_LOOKALIKES
                    )
                )
            ]
            if (
                len(found.get("EMAIL_ADDRESS", set())) < MAX_LOCATIONS_PER_CODE
                and public_email_matches
            ):
                found["EMAIL_ADDRESS"].add(location)

            lower = segment.lower()
            if (
                len(found.get("POTENTIALLY_PRIVATE_URL", set())) < MAX_LOCATIONS_PER_CODE
                and URL.search(segment)
                and any(
                    word in lower
                    for word in ("reviewer", "private", "access_token")
                )
            ):
                found["POTENTIALLY_PRIVATE_URL"].add(location)
    return found


def merge_locations(
    target: dict[str, set[str]], source: dict[str, set[str]]
) -> None:
    for code, locations in source.items():
        target[code].update(locations)


def check_git_state(
    root: Path, git_executable: str, reporter: Reporter
) -> None:
    status = run_git(root, git_executable, ["status", "--porcelain=v1"])
    dirty_lines = [line for line in status.stdout.splitlines() if line]
    if dirty_lines:
        reporter.add(
            "WARN",
            "WORKTREE_DIRTY",
            f"The working tree has {len(dirty_lines)} uncommitted change(s); commit or intentionally discard them before publication.",
        )
    else:
        reporter.add("PASS", "WORKTREE_CLEAN", "The working tree is clean.")

    upstream = run_git(
        root,
        git_executable,
        ["rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{upstream}"],
        check=False,
    )
    if upstream.returncode != 0:
        reporter.add("WARN", "NO_UPSTREAM", "The current branch has no upstream branch.")
        return
    counts = run_git(
        root,
        git_executable,
        ["rev-list", "--left-right", "--count", "@{upstream}...HEAD"],
    ).stdout.decode("ascii", errors="replace").strip().split()
    behind, ahead = (int(counts[0]), int(counts[1]))
    if behind or ahead:
        reporter.add(
            "WARN",
            "UPSTREAM_DIVERGENCE",
            f"Local branch is {ahead} commit(s) ahead and {behind} commit(s) behind its recorded upstream; fetch/push as appropriate.",
        )
    else:
        reporter.add("PASS", "UPSTREAM_SYNC", "HEAD matches the recorded upstream branch.")


def tracked_files(root: Path, git_executable: str) -> list[str]:
    result = run_git(root, git_executable, ["ls-files", "-z"])
    return decode_paths(result.stdout)


def check_community_files(root: Path, reporter: Reporter) -> None:
    if (root / "README.md").is_file():
        reporter.add("PASS", "README_PRESENT", "README.md is present.")
    else:
        reporter.add("BLOCKER", "README_MISSING", "Add a root README before publication.")

    license_files = [
        path for path in root.iterdir()
        if path.is_file() and path.name.lower().startswith(("license", "licence"))
    ]
    if license_files:
        reporter.add("PASS", "LICENSE_PRESENT", "A root license file is present.")
    else:
        reporter.add(
            "BLOCKER",
            "LICENSE_MISSING",
            "Choose and add a license before describing the repository as reusable open-source software.",
        )

    if (root / "CITATION.cff").is_file():
        reporter.add("PASS", "CITATION_PRESENT", "CITATION.cff is present.")
    else:
        reporter.add(
            "WARN", "CITATION_MISSING",
            "Add CITATION.cff so GitHub can display a machine-readable citation."
        )

    if (root / "SECURITY.md").is_file() or (root / ".github" / "SECURITY.md").is_file():
        reporter.add("PASS", "SECURITY_POLICY_PRESENT", "A security policy is present.")
    else:
        reporter.add(
            "WARN", "SECURITY_POLICY_MISSING",
            "Add SECURITY.md or enable an appropriate private vulnerability-reporting route."
        )

    readme_path = root / "README.md"
    if readme_path.is_file():
        readme = readme_path.read_text(encoding="utf-8-sig", errors="replace")
        required_references = {
            "PUBLIC_REPOSITORY_URL": "github.com/bomacchih/maize_shoot_data_process_v2",
            "ZENODO_RECORD": "zenodo.org/records/22058284",
            "PUBLICATION_DOI": "10.1111/pbi.70515",
        }
        missing = [code for code, needle in required_references.items() if needle not in readme]
        if missing:
            reporter.add(
                "WARN", "README_RELEASE_REFERENCES",
                "README is missing one or more repository/data/publication references: " + ", ".join(missing)
            )
        else:
            reporter.add(
                "PASS", "README_RELEASE_REFERENCES",
                "README records the public repository URL, Zenodo record, and publication DOI."
            )


def check_ignore_rules(root: Path, reporter: Reporter) -> None:
    path = root / ".gitignore"
    if not path.is_file():
        reporter.add("BLOCKER", "GITIGNORE_MISSING", "Add .gitignore before publication.")
        return
    lines = {
        line.strip().lower()
        for line in path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    required_groups = {
        "RDS": {"*.rds"},
        "HDF5": {"*.h5", "*.h5ad"},
        "LOOM": {"*.loom"},
        "FASTQ": {"*.fastq", "*.fastq.gz", "*.fq", "*.fq.gz"},
        "CLOUPE": {"*.cloupe"},
        "R_HISTORY": {".rhistory"},
        "LOCAL_PATHS": {"config/local_paths.yaml"},
    }
    missing = [
        label for label, alternatives in required_groups.items()
        if not alternatives.intersection(lines)
    ]
    if missing:
        reporter.add(
            "BLOCKER", "GITIGNORE_GAPS",
            "Missing ignore protection for: " + ", ".join(missing)
        )
    else:
        reporter.add(
            "PASS", "GITIGNORE_SAFEGUARDS",
            "Large data, R history, and machine-specific paths are ignored."
        )


def check_current_files(
    root: Path,
    paths: list[str],
    reporter: Reporter,
) -> None:
    sensitive = sorted(path for path in paths if has_sensitive_suffix(path))
    if sensitive:
        reporter.add(
            "BLOCKER", "TRACKED_SENSITIVE_FILES",
            "Sensitive credential or large analysis-data formats are tracked.", sensitive
        )
    else:
        reporter.add(
            "PASS", "NO_TRACKED_SENSITIVE_FILES",
            "No blocked credential or large analysis-data formats are tracked."
        )

    large_warnings: list[str] = []
    large_blockers: list[str] = []
    for relative in paths:
        path = root / relative
        if not path.is_file():
            continue
        size = path.stat().st_size
        label = f"{relative} ({size / MIB:.1f} MiB)"
        if size >= 100 * MIB:
            large_blockers.append(label)
        elif size >= 10 * MIB:
            large_warnings.append(label)
    if large_blockers:
        reporter.add(
            "BLOCKER", "TRACKED_FILE_OVER_100_MIB",
            "GitHub blocks normal Git files of 100 MiB or larger.", large_blockers
        )
    if large_warnings:
        reporter.add(
            "WARN", "LARGE_TRACKED_FILES",
            "Review tracked files of 10 MiB or larger; consider Zenodo or release assets when they are generated outputs.",
            large_warnings,
        )
    if not large_blockers and not large_warnings:
        reporter.add("PASS", "TRACKED_FILE_SIZES", "No tracked file is 10 MiB or larger.")

    all_findings: dict[str, set[str]] = defaultdict(set)
    for relative in paths:
        if relative == AUDIT_SCRIPT or not is_text_path(relative):
            continue
        path = root / relative
        if not path.is_file():
            continue
        try:
            text = path.read_text(encoding="utf-8-sig", errors="replace")
        except OSError:
            continue
        merge_locations(all_findings, classify_text_findings(text, relative))
    report_text_findings(all_findings, reporter, history=False)


def report_text_findings(
    findings: dict[str, set[str]],
    reporter: Reporter,
    *,
    history: bool,
) -> None:
    prefix = "HISTORY" if history else "CURRENT"
    secret_codes = set(SECRET_PATTERNS)
    secret_locations = sorted(
        location
        for code in secret_codes
        for location in findings.get(code, set())
    )
    if secret_locations:
        reporter.add(
            "BLOCKER", f"{prefix}_POTENTIAL_SECRETS",
            "Potential credential material was detected; values are redacted. Inspect and rotate/remove any genuine credential.",
            secret_locations,
        )
    else:
        reporter.add(
            "PASS", f"{prefix}_SECRET_SCAN",
            "No supported high-confidence credential pattern was detected."
        )

    private_urls = sorted(findings.get("POTENTIALLY_PRIVATE_URL", set()))
    if private_urls:
        reporter.add(
            "BLOCKER", f"{prefix}_POTENTIALLY_PRIVATE_URLS",
            "URLs near reviewer/private-access terms require manual review; URL values are redacted.",
            private_urls,
        )

    personal_paths = sorted(findings.get("PERSONAL_ABSOLUTE_PATH", set()))
    if personal_paths:
        reporter.add(
            "WARN", f"{prefix}_PERSONAL_PATHS",
            "Personal absolute filesystem paths are present; decide whether to replace them with variables or examples.",
            personal_paths,
        )
    else:
        reporter.add(
            "PASS", f"{prefix}_PERSONAL_PATHS",
            "No non-placeholder personal home-directory path was detected."
        )

    emails = sorted(findings.get("EMAIL_ADDRESS", set()))
    if emails:
        reporter.add(
            "WARN", f"{prefix}_EMAIL_ADDRESSES",
            "Email addresses are present; verify that each is intended for public release.",
            emails,
        )


def history_blob_inventory(
    root: Path,
    git_executable: str,
) -> tuple[list[tuple[str, str, int]], dict[str, str]]:
    objects = run_git(root, git_executable, ["rev-list", "--objects", "--all"])
    path_for_object: dict[str, str] = {}
    object_ids: list[str] = []
    for raw_line in objects.stdout.decode("utf-8", errors="surrogateescape").splitlines():
        parts = raw_line.split(" ", 1)
        object_id = parts[0]
        object_ids.append(object_id)
        if len(parts) == 2:
            path_for_object.setdefault(object_id, parts[1])
    payload = ("\n".join(object_ids) + "\n").encode("ascii")
    checked = run_git(
        root,
        git_executable,
        ["cat-file", "--batch-check=%(objectname) %(objecttype) %(objectsize)"],
        input_bytes=payload,
    )
    blobs: list[tuple[str, str, int]] = []
    for line in checked.stdout.decode("ascii", errors="replace").splitlines():
        object_id, object_type, size_text = line.split()
        if object_type == "blob":
            blobs.append((object_id, path_for_object.get(object_id, "<unknown>"), int(size_text)))
    return blobs, path_for_object


def scan_history(
    root: Path,
    git_executable: str,
    reporter: Reporter,
    maximum_text_blob_bytes: int,
) -> None:
    blobs, _ = history_blob_inventory(root, git_executable)
    sensitive_paths = sorted({path for _, path, _ in blobs if has_sensitive_suffix(path)})
    if sensitive_paths:
        reporter.add(
            "BLOCKER", "HISTORY_SENSITIVE_FILES",
            "Sensitive credential or large analysis-data formats occur in reachable Git history.",
            sensitive_paths,
        )
    else:
        reporter.add(
            "PASS", "HISTORY_SENSITIVE_FILES",
            "No blocked credential or large analysis-data format occurs in reachable history."
        )

    history_large_warnings: list[str] = []
    history_large_blockers: list[str] = []
    for _, path, size in blobs:
        label = f"{path} ({size / MIB:.1f} MiB)"
        if size >= 100 * MIB:
            history_large_blockers.append(label)
        elif size >= 10 * MIB:
            history_large_warnings.append(label)
    if history_large_blockers:
        reporter.add(
            "BLOCKER", "HISTORY_BLOB_OVER_100_MIB",
            "Reachable Git history contains blobs of 100 MiB or larger.",
            history_large_blockers,
        )
    if history_large_warnings:
        reporter.add(
            "WARN", "HISTORY_LARGE_BLOBS",
            "Reachable history contains blobs of 10 MiB or larger.",
            history_large_warnings,
        )
    if not history_large_blockers and not history_large_warnings:
        reporter.add("PASS", "HISTORY_BLOB_SIZES", "No reachable blob is 10 MiB or larger.")

    candidates = [
        (object_id, path, size)
        for object_id, path, size in blobs
        if path != AUDIT_SCRIPT and is_text_path(path) and size <= maximum_text_blob_bytes
    ]
    findings: dict[str, set[str]] = defaultdict(set)
    command = [
        git_executable,
        "-c",
        f"safe.directory={root.as_posix()}",
        "cat-file",
        "--batch",
    ]
    process = subprocess.Popen(
        command,
        cwd=root,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    assert process.stdin is not None and process.stdout is not None
    try:
        for object_id, path, _ in candidates:
            process.stdin.write((object_id + "\n").encode("ascii"))
            process.stdin.flush()
            header = process.stdout.readline().decode("ascii", errors="replace").strip()
            header_parts = header.split()
            if len(header_parts) != 3 or header_parts[1] != "blob":
                raise AuditError(f"Unexpected git cat-file response for {object_id[:12]}")
            size = int(header_parts[2])
            content = process.stdout.read(size)
            process.stdout.read(1)
            text = content.decode("utf-8", errors="replace")
            prefix = f"{path}@blob-{object_id[:12]}"
            merge_locations(findings, classify_text_findings(text, prefix))
    finally:
        if process.stdin:
            process.stdin.close()
        process.wait(timeout=30)
    if process.returncode != 0:
        raise AuditError("git cat-file history scan failed")
    report_text_findings(findings, reporter, history=True)


def write_json_report(path: Path, root: Path, reporter: Reporter) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "project_root": ".",
        "summary": {
            "passed": reporter.count("PASS"),
            "warnings": reporter.count("WARN"),
            "blockers": reporter.count("BLOCKER"),
        },
        "findings": [asdict(finding) for finding in reporter.findings],
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def write_markdown_report(path: Path, root: Path, reporter: Reporter) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    status = "READY" if reporter.count("BLOCKER") == 0 else "NOT READY"
    lines = [
        "# Public-readiness audit",
        "",
        f"- Status: **{status}**",
        f"- Repository: `{root.name}`",
        f"- Generated: {datetime.now(timezone.utc).isoformat()}",
        f"- Summary: {reporter.count('PASS')} passed, {reporter.count('WARN')} warnings, {reporter.count('BLOCKER')} blockers",
        "",
        "Potential secret values are redacted; locations identify where manual review is required.",
        "",
    ]
    for severity in ("BLOCKER", "WARN", "PASS"):
        lines.extend([f"## {severity.title()}", ""])
        matching = [finding for finding in reporter.findings if finding.severity == severity]
        if not matching:
            lines.extend(["None.", ""])
            continue
        for finding in matching:
            lines.append(f"### {finding.code}")
            lines.extend(["", finding.message, ""])
            for location in finding.locations:
                lines.append(f"- `{location}`")
            if finding.locations:
                lines.append("")
    path.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    reporter = Reporter()
    try:
        top_level = run_git(root, args.git_executable, ["rev-parse", "--show-toplevel"])
        root = Path(top_level.stdout.decode("utf-8", errors="replace").strip()).resolve()
        paths = tracked_files(root, args.git_executable)
        check_git_state(root, args.git_executable, reporter)
        check_community_files(root, reporter)
        check_ignore_rules(root, reporter)
        check_current_files(root, paths, reporter)
        if args.skip_history:
            reporter.add("WARN", "HISTORY_SCAN_SKIPPED", "Full-history scanning was skipped.")
        else:
            scan_history(
                root,
                args.git_executable,
                reporter,
                maximum_text_blob_bytes=max(1, int(args.history_max_blob_mb * MIB)),
            )
    except (AuditError, OSError, subprocess.SubprocessError) as exc:
        reporter.add("BLOCKER", "AUDIT_INCOMPLETE", str(exc))

    reporter.print_console()
    if args.json_report:
        write_json_report(args.json_report.resolve(), root, reporter)
        print(f"JSON report: {args.json_report.resolve()}")
    if args.markdown_report:
        write_markdown_report(args.markdown_report.resolve(), root, reporter)
        print(f"Markdown report: {args.markdown_report.resolve()}")

    if reporter.count("BLOCKER"):
        return 1
    if args.strict and reporter.count("WARN"):
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
