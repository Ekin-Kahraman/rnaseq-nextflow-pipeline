"""Check paired-end fixture semantics independently against the reference sequence."""

import gzip
import random
import tempfile
import unittest
from pathlib import Path

from create_test_data import write_fastq_pair


class PairedFixtureTests(unittest.TestCase):
    def test_mates_share_an_inward_facing_reference_fragment(self):
        rng = random.Random(2026)
        reference = "".join(rng.choice("ACGT") for _ in range(30000))
        with tempfile.TemporaryDirectory() as directory:
            write_fastq_pair("paired", 500, directory, reference)
            records = []
            for mate in [1, 2]:
                with gzip.open(Path(directory) / f"paired_{mate}.fastq.gz", "rt") as handle:
                    lines = handle.read().splitlines()
                self.assertEqual(len(lines), 500 * 4)
                records.append([lines[i:i + 4] for i in range(0, len(lines), 4)])
            complement = str.maketrans("ACGT", "TGCA")
            for i, (r1, r2) in enumerate(zip(*records)):
                self.assertEqual(r1[0], f"@paired.{i}/R1")
                self.assertEqual(r2[0], f"@paired.{i}/R2")
                self.assertEqual(r1[2], "+")
                self.assertEqual(r2[2], "+")
                self.assertEqual(len(r1[1]), 100)
                self.assertEqual(len(r2[1]), 100)
                self.assertEqual(len(r1[3]), 100)
                self.assertEqual(len(r2[3]), 100)
                start = reference.find(r1[1])
                end_start = reference.find(r2[1].translate(complement)[::-1])
                self.assertGreaterEqual(start, 0)
                self.assertGreaterEqual(end_start, 0)
                self.assertEqual(reference.find(r1[1], start + 1), -1)
                self.assertEqual(reference.find(r2[1].translate(complement)[::-1], end_start + 1), -1)
                self.assertGreaterEqual(end_start - start, 150)
                self.assertLessEqual(end_start - start, 300)

    def test_gzip_output_is_reproducible_for_the_same_seed(self):
        reference = "ACGT" * 1000
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            state = random.getstate()
            try:
                random.seed(42)
                write_fastq_pair("paired", 20, directory, reference)
                first = [(path / f"paired_{mate}.fastq.gz").read_bytes() for mate in [1, 2]]
                random.seed(42)
                write_fastq_pair("paired", 20, directory, reference)
                second = [(path / f"paired_{mate}.fastq.gz").read_bytes() for mate in [1, 2]]
            finally:
                random.setstate(state)
            self.assertEqual(first, second)


if __name__ == "__main__":
    unittest.main()
