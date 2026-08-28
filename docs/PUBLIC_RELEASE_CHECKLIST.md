# Public release checklist

Run the automated audit from the repository root before changing GitHub
visibility:

```bash
python scripts/python/validate_public_readiness.py \
  --json-report results/logs/public_readiness_report.json \
  --markdown-report results/logs/public_readiness_report.md
```

If Git is not on `PATH`, supply its executable with `--git-executable` or set
the `GIT_EXE` environment variable.

The audit reports locations but never prints matching credential values. It
checks the current tracked tree and all reachable Git history for:

- potential credentials and private/reviewer URLs;
- personal absolute filesystem paths and public email addresses;
- tracked or historical RDS, HDF5, loom, BAM, FASTQ, `.cloupe`, key, and
  environment files;
- GitHub's 100 MiB file limit and files over 10 MiB that merit review;
- `.gitignore` safeguards for large or machine-specific data;
- repository/upstream state; and
- README, license, citation, and security-policy readiness.

## Manual checks

Automation cannot decide the following items:

- Select a software/content license that all rights holders approve.
- Confirm that every author and contributor approves the public release.
- Confirm that manuscript, reviewer, embargo, human-subject, and third-party
  data restrictions permit publication.
- Review every email address and personal filesystem path reported by the
  audit, including historical occurrences.
- Confirm that large generated reports and figures belong in Git rather than
  Zenodo or GitHub release assets.
- Verify the Zenodo record, DOI, accessions, and manuscript citation from a
  clean clone.
- Run the configuration, R, Python, and Seurat-object validation suites.
- Inspect GitHub Actions history and logs before changing visibility; those
  logs become publicly visible with the repository.
- Enable appropriate GitHub secret scanning, push protection, and private
  vulnerability reporting after publication.

## Release gate

Do not make the repository public while the audit reports a `BLOCKER`.
Warnings require a documented decision but do not necessarily prevent release.
After addressing findings, commit the changes, push them, rerun the audit from
a clean clone, and only then change repository visibility.

