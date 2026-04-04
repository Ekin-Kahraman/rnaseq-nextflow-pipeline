#!/usr/bin/env python3
"""Generate minimal synthetic test data for pipeline verification.

Creates tiny FASTQ files, a minimal genome, and a gene annotation that
are just large enough for each pipeline step to run without error.
"""

import gzip
import random
import subprocess
import os

random.seed(42)

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
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
    """Write paired FASTQ files with reads sampled from the genome."""
    max_start = len(genome_seq) - READ_LEN * 2 - 200
    if max_start < 0:
        max_start = 0
    for read_num in [1, 2]:
        path = os.path.join(outdir, f"{sample_id}_{read_num}.fastq.gz")
        with gzip.open(path, "wt") as f:
            for i in range(n_reads):
                # Sample a position from the genome
                start = random.randint(0, max_start)
                if read_num == 1:
                    seq = genome_seq[start : start + READ_LEN]
                else:
                    # Mate reads from ~200bp downstream
                    mate_start = start + READ_LEN + random.randint(50, 200)
                    mate_start = min(mate_start, len(genome_seq) - READ_LEN)
                    seq = genome_seq[mate_start : mate_start + READ_LEN]
                    # Reverse complement for paired-end
                    comp = str.maketrans("ACGT", "TGCA")
                    seq = seq.translate(comp)[::-1]
                f.write(f"@{sample_id}.{i}/R{read_num}\n")
                f.write(seq + "\n")
                f.write("+\n")
                f.write(QUAL + "\n")
    print(f"  {sample_id}: {n_reads} read pairs (sampled from genome)")


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
        f.write("sample_id,fastq_1,fastq_2,condition\n")
        for sample_id, condition in samples:
            f.write(
                f"{sample_id},"
                f"{outdir}/{sample_id}_1.fastq.gz,"
                f"{outdir}/{sample_id}_2.fastq.gz,"
                f"{condition}\n"
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
        ("test_pos1", "positive"),
        ("test_pos2", "positive"),
        ("test_neg1", "negative"),
        ("test_neg2", "negative"),
    ]

    genome_path = write_genome(TEST_DIR)

    # Read back the genome for read sampling
    with open(genome_path) as f:
        genome_seq = "".join(line.strip() for line in f if not line.startswith(">"))

    for sample_id, _ in samples:
        write_fastq_pair(sample_id, N_READS, TEST_DIR, genome_seq)
    write_gtf(TEST_DIR)
    write_samplesheet(samples, TEST_DIR)
    build_hisat2_index(genome_path, TEST_DIR)

    print("Done.")


if __name__ == "__main__":
    main()
