suppressPackageStartupMessages(library(DESeq2))
source("scripts/deseq2_fallback.R")

counts <- matrix(
  c(10L, 12L, 8L, 11L, 20L, 18L, 22L, 19L),
  nrow = 2,
  dimnames = list(c("gene_a", "gene_b"), c("a1", "a2", "b1", "b2"))
)
sample_data <- S4Vectors::DataFrame(
  condition = factor(c("control", "control", "treated", "treated"))
)
dds <- DESeqDataSetFromMatrix(counts, sample_data, design = ~ condition)

forced_failure <- function(dds) {
  stop("intentional test failure")
}

error <- tryCatch(
  run_deseq_with_test_fallback(dds, allow_test_fallback = FALSE, de_function = forced_failure),
  error = identity
)

stopifnot(inherits(error, "error"))
stopifnot(grepl("only in the synthetic test profile", conditionMessage(error), fixed = TRUE))
stopifnot(grepl("intentional test failure", conditionMessage(error), fixed = TRUE))
message("Passed: real-profile DESeq failures are reported and do not trigger a model fallback.")
