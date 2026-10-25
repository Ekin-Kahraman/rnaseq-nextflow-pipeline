run_deseq_with_test_fallback <- function(dds, allow_test_fallback = FALSE, de_function = DESeq) {
  tryCatch(
    de_function(dds),
    error = function(error) {
      if (!isTRUE(allow_test_fallback)) {
        stop(
          paste(
            "DESeq failed. Gene-wise dispersion fallback is enabled only in the synthetic test profile.",
            "Original error:", conditionMessage(error)
          ),
          call. = FALSE
        )
      }

      message(
        "Synthetic test profile only: using gene-wise dispersion estimates after DESeq failed. ",
        "Reason: ", conditionMessage(error)
      )
      dds <- estimateSizeFactors(dds, type = "poscounts")
      dds <- estimateDispersionsGeneEst(dds)
      dispersions(dds) <- mcols(dds)$dispGeneEst
      nbinomWaldTest(dds)
    }
  )
}
