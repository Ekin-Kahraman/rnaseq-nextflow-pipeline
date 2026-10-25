#!/usr/bin/env python3
"""Validate the minimal artefacts expected from a completed pipeline run."""

from __future__ import annotations

import argparse
import csv
import math
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
DESEQ_NUMERIC_COLUMNS = (
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj",
)


def require_file(path: Path) -> None:
    if not path.exists():
        raise SystemExit(f"Missing expected output: {path}")
    if path.stat().st_size == 0:
        raise SystemExit(f"Output is empty: {path}")


def read_sample_ids(path: Path) -> set[str]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"sample_id", "fastq_1", "fastq_2", "condition"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise SystemExit(f"Sample sheet must include columns {sorted(required)}: {path}")
        rows = list(reader)
        sample_ids = [(row.get("sample_id") or "").strip() for row in rows]

    if not sample_ids or any(not sample_id for sample_id in sample_ids):
        raise SystemExit(f"Sample sheet contains an empty or missing sample ID: {path}")
    if len(sample_ids) != len(set(sample_ids)):
        raise SystemExit(f"Sample sheet contains duplicate sample IDs: {path}")
    return set(sample_ids)


def validate_counts(path: Path, expected_samples: set[str] | None = None) -> set[str]:
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        rows = [row for row in reader if row and not row[0].startswith("#")]

    if len(rows) < 2:
        raise SystemExit(f"Count matrix has no gene rows: {path}")

    header = rows[0]
    if len(header) < 7 or header[0] != "Geneid":
        raise SystemExit(f"Count matrix header is not featureCounts format: {path}")

    sample_ids = [Path(column.replace("\\", "/")).name.removesuffix(".bam") for column in header[6:]]
    if any(not sample_id for sample_id in sample_ids) or len(sample_ids) != len(set(sample_ids)):
        raise SystemExit(f"Count matrix has empty or duplicate sample columns: {path}")
    if expected_samples is not None and set(sample_ids) != expected_samples:
        missing = sorted(expected_samples.difference(sample_ids))
        extra = sorted(set(sample_ids).difference(expected_samples))
        raise SystemExit(f"Count matrix samples do not match the sample sheet (missing={missing}, extra={extra}): {path}")

    non_zero_rows = 0
    gene_ids: set[str] = set()
    for fields in rows[1:]:
        if len(fields) != len(header):
            raise SystemExit(f"Malformed count row in {path}: {' '.join(fields)[:120]}")
        gene_id = fields[0].strip()
        if not gene_id or gene_id in gene_ids:
            raise SystemExit(f"Count matrix has an empty or duplicate gene ID: {gene_id!r} in {path}")
        gene_ids.add(gene_id)
        try:
            counts = [int(value) for value in fields[6:]]
        except ValueError as error:
            raise SystemExit(f"Count matrix contains a non-integer count in {path}: {error}") from error
        if any(value < 0 for value in counts):
            raise SystemExit(f"Count matrix contains a negative count in {path}: {gene_id}")
        if any(value > 0 for value in counts):
            non_zero_rows += 1

    if non_zero_rows == 0:
        raise SystemExit(f"Count matrix contains no non-zero gene counts: {path}")
    return gene_ids


def validate_deseq(path: Path, count_gene_ids: set[str] | None = None) -> set[str]:
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
    if len(genes) != len(rows):
        raise SystemExit(f"DESeq2 results contain empty or duplicate gene IDs: {path}")
    if count_gene_ids is not None and not genes.issubset(count_gene_ids):
        unknown = sorted(genes.difference(count_gene_ids))[:10]
        raise SystemExit(f"DESeq2 results contain genes absent from the count matrix: {unknown}")

    for row in rows:
        for column in DESEQ_NUMERIC_COLUMNS:
            value = row[column].strip()
            if value in {"", "NA", "NaN"}:
                continue
            try:
                number = float(value)
            except ValueError as error:
                raise SystemExit(f"DESeq2 column {column} contains a non-numeric value: {value!r}") from error
            if not math.isfinite(number):
                raise SystemExit(f"DESeq2 column {column} contains a non-finite value: {value!r}")
            if column in {"pvalue", "padj"} and not 0 <= number <= 1:
                raise SystemExit(f"DESeq2 column {column} is outside [0, 1]: {value!r}")
    return genes


def validate_multiqc(path: Path) -> None:
    text = path.read_text(errors="ignore").lower()
    if "<html" not in text or "multiqc" not in text:
        raise SystemExit(f"MultiQC report does not look like an HTML report: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("outdir", type=Path, help="Pipeline output directory")
    parser.add_argument(
        "--samplesheet",
        type=Path,
        help="Optional sample sheet to verify that the count matrix contains the same sample IDs.",
    )
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

    expected_samples = read_sample_ids(args.samplesheet) if args.samplesheet else None
    count_gene_ids = validate_counts(args.outdir / "counts/gene_counts.txt", expected_samples)
    validate_deseq(args.outdir / "deseq2/deseq2_results.csv", count_gene_ids)
    validate_multiqc(args.outdir / "multiqc/multiqc_report.html")

    print(f"Validated RNA-seq pipeline outputs in {args.outdir}")


if __name__ == "__main__":
    main()
