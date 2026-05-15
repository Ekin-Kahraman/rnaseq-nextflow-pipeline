# Security policy

## Supported versions

Security fixes are made on the default branch. If releases are cut in future, only the latest release line will be supported.

## Reporting a vulnerability

This pipeline processes user-supplied FASTQ, CSV, GTF, reference-index and path inputs. If you find a vulnerability such as unsafe path handling, unintended file overwrite, container escape risk, or credential leakage in cloud execution, do not open a public issue.

Email: evk23umu@uea.ac.uk
Subject: `rnaseq-nextflow security: <short description>`

Please include:
- A minimal reproduction command or config
- Operating system and Nextflow version
- Whether Docker, Singularity, local, or AWS Batch was used
- Expected impact

I aim to respond within 7 days. Public correctness bugs, failed runs, and feature requests should be opened as GitHub issues.

## Data handling

Do not attach patient-level FASTQ files, clinical metadata, access keys, S3 bucket names, or private sample sheets to public issues. Use synthetic data or redacted paths wherever possible.
