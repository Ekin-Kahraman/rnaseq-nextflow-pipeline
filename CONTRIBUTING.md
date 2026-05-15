# Contributing

Contributions are welcome when they improve reproducibility, correctness, portability, or documentation.

## Development setup

```bash
git clone https://github.com/Ekin-Kahraman/rnaseq-nextflow-pipeline.git
cd rnaseq-nextflow-pipeline
python test/create_test_data.py
nextflow run main.nf -profile test,docker --outdir test_results
python scripts/validate_outputs.py test_results
```

Use `-profile test` for pull requests. It builds a tiny synthetic reference and verifies the full software path without downloading controlled or patient-derived data.

## Pull request checklist

- The test profile completes locally.
- `python scripts/validate_outputs.py <outdir>` passes.
- New parameters are documented in `README.md` and `nextflow.config`.
- Changes affecting cloud execution update `docs/cloud.md`.
- Output layout changes are reflected in the README output tree.

## Scope

Preferred contributions:
- Additional aligner or quantifier options with tests
- More robust schema validation
- Cloud execution improvements
- Documentation for real datasets
- Bug fixes in output validation

Out of scope:
- Committing real FASTQ files, large references, or private metadata
- Unpinned custom containers without a clear reason
- Changes that silently alter biological interpretation without documenting the effect

## Licence

By contributing, you agree that your contributions are licensed under the MIT licence.
