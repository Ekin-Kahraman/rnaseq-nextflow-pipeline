#!/usr/bin/env python3
"""Validate the minimal artefacts expected from a completed pipeline run."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


REQUIRED_FILES = [
    "counts/gene_counts.txt",
    "counts/gene_counts.txt.summary",
    "deseq2/deseq2_results.csv",
    "deseq2/volcano_plot.png",
    "deseq2/pca_plot.png",
    "multiqc/multiqc_report.html",
]

RUN_METADATA_FILES = [
    "pipeline_info/report.html",
    "pipeline_info/timeline.html",
    "pipeline_info/trace.txt",
    "pipeline_info/dag.dot",
]

DESEQ_COLUMNS = {
    "gene",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
}


def require_file(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"Missing expected output: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"Output is empty: {path}")


def validate_counts(path: Path) -> None:
    lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
    if len(lines) < 2:
        raise SystemExit(f"Count matrix has no gene rows: {path}")

    header = lines[0].split("\t")
    if len(header) < 7 or header[0] != "Geneid":
        raise SystemExit(f"Count matrix header is not featureCounts format: {path}")

    non_zero_rows = 0
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != len(header):
            raise SystemExit(f"Malformed count row in {path}: {line[:120]}")
        counts = [int(value) for value in fields[6:]]
        if any(value > 0 for value in counts):
            non_zero_rows += 1

    if non_zero_rows == 0:
        raise SystemExit(f"Count matrix contains no non-zero gene counts: {path}")


def validate_deseq(path: Path) -> None:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise SystemExit(f"DESeq2 results have no header: {path}")

        missing = DESEQ_COLUMNS.difference(reader.fieldnames)
        if missing:
            raise SystemExit(f"DESeq2 results missing columns {sorted(missing)}: {path}")

        rows = list(reader)
    if not rows:
        raise SystemExit(f"DESeq2 results have no gene rows: {path}")

    genes = {row["gene"] for row in rows if row.get("gene")}
    if not genes:
        raise SystemExit(f"DESeq2 results have no named genes: {path}")


def validate_multiqc(path: Path) -> None:
    text = path.read_text(errors="ignore").lower()
    if "<html" not in text or "multiqc" not in text:
        raise SystemExit(f"MultiQC report does not look like an HTML report: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Pipeline output directory")
    parser.add_argument(
        "--skip-run-metadata",
        action="store_true",
        help="Do not require Nextflow report, timeline, trace and DAG files",
    )
    args = parser.parse_args()

    required = list(REQUIRED_FILES)
    if not args.skip_run_metadata:
        required.extend(RUN_METADATA_FILES)

    for relpath in required:
        require_file(args.outdir / relpath)

    validate_counts(args.outdir / "counts/gene_counts.txt")
    validate_deseq(args.outdir / "deseq2/deseq2_results.csv")
    validate_multiqc(args.outdir / "multiqc/multiqc_report.html")

    print(f"Validated RNA-seq pipeline outputs in {args.outdir}")


if __name__ == "__main__":
    main()
