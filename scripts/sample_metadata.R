read_sample_metadata <- function(path, sample_ids) {
  # Identifiers are literal strings: 01, 1 and NA can identify different people.
  samples <- read.csv(
    path, colClasses = "character", na.strings = character(),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  if (!all(c("sample_id", "condition") %in% names(samples))) {
    stop("Sample metadata requires sample_id and condition columns.")
  }
  if (anyDuplicated(samples$sample_id) || any(!nzchar(samples$sample_id))) {
    stop("Sample metadata requires unique, non-empty sample IDs.")
  }
  if (anyDuplicated(sample_ids) || !setequal(sample_ids, samples$sample_id)) {
    stop("Count columns and sample metadata must contain the same unique sample IDs.")
  }
  samples <- samples[match(sample_ids, samples$sample_id), , drop = FALSE]
  rownames(samples) <- sample_ids
  samples
}
