#!/usr/bin/env python3
"""Generate minimal synthetic test data for pipeline verification.

Creates tiny FASTQ files, a minimal genome, and a gene annotation that
are just large enough for each pipeline step to run without error.
"""

from contextlib import ExitStack
import gzip
import io
import random
import subprocess
import os
from pathlib import Path

random.seed(42)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = os.path.dirname(TEST_DIR)
N_READS = 2000
READ_LEN = 100
GENE_LEN = 500
N_GENES = 50
SPACER = 100

BASES = "ACGT"
QUAL = "I" * READ_LEN  # Phred 40


def random_seq(length):
    return "".join(random.choice(BASES) for _ in range(length))


def write_fastq_pair(sample_id, n_reads, outdir, genome_seq):
    """Write inward-facing mates from the same 250-400 bp reference fragment."""
    min_fragment = READ_LEN * 2 + 50
    max_fragment = READ_LEN * 2 + 200
    if len(genome_seq) < max_fragment:
        raise ValueError(f"Synthetic genome must be at least {max_fragment} bp long.")
    complement = str.maketrans("ACGT", "TGCA")

    with ExitStack() as stack:
        mates = []
        for read_num in [1, 2]:
            path = os.path.join(outdir, f"{sample_id}_{read_num}.fastq.gz")
            raw = stack.enter_context(open(path, "wb"))
            gz = gzip.GzipFile(
                filename=os.path.basename(path).removesuffix(".gz"),
                mode="wb", fileobj=raw, mtime=0,
            )
            mates.append(stack.enter_context(io.TextIOWrapper(gz, encoding="utf-8", newline="")))

        for i in range(n_reads):
            fragment_length = random.randint(min_fragment, max_fragment)
            start = random.randint(0, len(genome_seq) - fragment_length)
            # Both mates come from one fragment; R2 faces back towards R1.
            r1 = genome_seq[start : start + READ_LEN]
            mate_start = start + fragment_length - READ_LEN
            r2 = genome_seq[mate_start : mate_start + READ_LEN].translate(complement)[::-1]
            for read_num, (handle, seq) in enumerate(zip(mates, [r1, r2]), start=1):
                handle.write(f"@{sample_id}.{i}/R{read_num}\n{seq}\n+\n{QUAL}\n")
    print(f"  {sample_id}: {n_reads} read pairs (sampled from shared reference fragments)")


def write_genome(outdir):
    """Write a FASTA genome with enough space for multiple genes."""
    genome_len = N_GENES * (GENE_LEN + SPACER)
    path = os.path.join(outdir, "genome.fa")
    seq = random_seq(genome_len)
    with open(path, "w") as f:
        f.write(">chr1\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i : i + 80] + "\n")
    print(f"  Genome: {genome_len} bp")
    return path


def write_gtf(outdir):
    """Write GTF annotation with multiple genes."""
    path = os.path.join(outdir, "genes.gtf")
    with open(path, "w") as f:
        for i in range(N_GENES):
            name = f"GENE{i+1:03d}"
            start = i * (GENE_LEN + SPACER) + 1
            end = start + GENE_LEN - 1
            f.write(
                f'chr1\ttest\tgene\t{start}\t{end}\t.\t+\t.\t'
                f'gene_id "{name}"; gene_name "{name}";\n'
            )
            f.write(
                f'chr1\ttest\texon\t{start}\t{end}\t.\t+\t.\t'
                f'gene_id "{name}"; gene_name "{name}";\n'
            )
    print(f"  GTF: {N_GENES} genes")
    return path


def write_samplesheet(samples, outdir):
    """Write test samplesheet."""
    path = os.path.join(outdir, "samplesheet.csv")
    with open(path, "w") as f:
        f.write("sample_id,fastq_1,fastq_2,condition,donor\n")
        for sample_id, condition, donor in samples:
            fastq_1 = Path(outdir, f"{sample_id}_1.fastq.gz").relative_to(PROJECT_DIR)
            fastq_2 = Path(outdir, f"{sample_id}_2.fastq.gz").relative_to(PROJECT_DIR)
            f.write(
                f"{sample_id},"
                f"{fastq_1},"
                f"{fastq_2},"
                f"{condition},"
                f"{donor}\n"
            )
    print(f"  Samplesheet: {len(samples)} samples")


def build_hisat2_index(genome_path, outdir):
    """Build HISAT2 index from the test genome."""
    prefix = os.path.join(outdir, "genome")
    try:
        subprocess.run(
            ["hisat2-build", "-q", genome_path, prefix],
            check=True, capture_output=True,
        )
        print(f"  HISAT2 index built: {prefix}")
    except FileNotFoundError:
        print("  Warning: hisat2-build not found locally. Index must be built in Docker.")


def main():
    print("Generating test data...")

    samples = [
        ("test_pos1", "positive", "donor_1"),
        ("test_pos2", "positive", "donor_2"),
        ("test_neg1", "negative", "donor_1"),
        ("test_neg2", "negative", "donor_2"),
    ]

    genome_path = write_genome(TEST_DIR)

    # Read back the genome for read sampling
    with open(genome_path) as f:
        genome_seq = "".join(line.strip() for line in f if not line.startswith(">"))

    for sample_id, _, _ in samples:
        write_fastq_pair(sample_id, N_READS, TEST_DIR, genome_seq)
    write_gtf(TEST_DIR)
    write_samplesheet(samples, TEST_DIR)
    build_hisat2_index(genome_path, TEST_DIR)

    print("Done.")


if __name__ == "__main__":
    main()
