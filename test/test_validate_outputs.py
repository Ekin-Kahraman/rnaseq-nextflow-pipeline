import csv
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.validate_outputs import read_sample_ids, validate_counts, validate_deseq


class OutputValidationTests(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.root = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def write_counts(self, body: str) -> Path:
        path = self.root / "gene_counts.txt"
        path.write_text(
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/work/test_a.bam\t/work/test_b.bam\n"
            + body,
            encoding="utf-8",
        )
        return path

    def write_deseq(self, body: str) -> Path:
        path = self.root / "deseq2_results.csv"
        path.write_text(
            "gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj\n" + body,
            encoding="utf-8",
        )
        return path

    def test_counts_validate_sample_ids_and_gene_counts(self):
        path = self.write_counts("GENE1\tchr1\t1\t10\t+\t10\t5\t0\n")
        gene_ids = validate_counts(path, {"test_a", "test_b"})
        self.assertEqual(gene_ids, {"GENE1"})

    def test_counts_reject_mismatched_samples(self):
        path = self.write_counts("GENE1\tchr1\t1\t10\t+\t10\t5\t0\n")
        with self.assertRaisesRegex(SystemExit, "do not match"):
            validate_counts(path, {"test_a", "test_c"})

    def test_counts_reject_duplicate_sample_names(self):
        path = self.root / "gene_counts.txt"
        path.write_text(
            "Geneid\tChr\tStart\tEnd\tStrand\tLength\t/work/test_a.bam\t/test_a.bam\n"
            "GENE1\tchr1\t1\t10\t+\t10\t5\t0\n",
            encoding="utf-8",
        )
        with self.assertRaisesRegex(SystemExit, "duplicate sample columns"):
            validate_counts(path)

    def test_counts_reject_negative_and_non_integer_values(self):
        negative = self.write_counts("GENE1\tchr1\t1\t10\t+\t10\t-1\t0\n")
        with self.assertRaisesRegex(SystemExit, "negative count"):
            validate_counts(negative)

        non_integer = self.write_counts("GENE1\tchr1\t1\t10\t+\t10\t1.5\t0\n")
        with self.assertRaisesRegex(SystemExit, "non-integer"):
            validate_counts(non_integer)

    def test_counts_reject_duplicate_gene_ids(self):
        path = self.write_counts(
            "GENE1\tchr1\t1\t10\t+\t10\t5\t0\n"
            "GENE1\tchr1\t11\t20\t+\t10\t3\t2\n"
        )
        with self.assertRaisesRegex(SystemExit, "duplicate gene ID"):
            validate_counts(path)

    def test_deseq_rejects_genes_absent_from_counts(self):
        path = self.write_deseq("GENE2,10,1.2,0.2,6,0.001,0.01\n")
        with self.assertRaisesRegex(SystemExit, "absent from the count matrix"):
            validate_deseq(path, {"GENE1"})

    def test_deseq_rejects_invalid_probability(self):
        path = self.write_deseq("GENE1,10,1.2,0.2,6,1.2,0.01\n")
        with self.assertRaisesRegex(SystemExit, r"outside \[0, 1\]"):
            validate_deseq(path, {"GENE1"})

    def test_sample_sheet_ids_are_read_and_unique(self):
        path = self.root / "samplesheet.csv"
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=["sample_id", "fastq_1", "fastq_2", "condition"])
            writer.writeheader()
            writer.writerow({"sample_id": "test_a", "fastq_1": "a_1.fq.gz", "fastq_2": "a_2.fq.gz", "condition": "a"})
            writer.writerow({"sample_id": "test_b", "fastq_1": "b_1.fq.gz", "fastq_2": "b_2.fq.gz", "condition": "b"})
        self.assertEqual(read_sample_ids(path), {"test_a", "test_b"})


if __name__ == "__main__":
    unittest.main()
