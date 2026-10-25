source("scripts/sample_metadata.R")

sheet <- tempfile(fileext = ".csv")
metadata <- data.frame(
  sample_id = c("001", "1", "002", "2", "NA", "003"),
  condition = rep(c("NA", "treated"), 3),
  donor = rep(c("01", "1", "2"), each = 2)
)
write.csv(metadata, sheet, row.names = FALSE, quote = FALSE)
ids <- rev(metadata$sample_id)
parsed <- read_sample_metadata(sheet, ids)
stopifnot(identical(parsed$sample_id, ids))
stopifnot(identical(parsed$donor, rev(metadata$donor)))
stopifnot(identical(parsed$condition, rev(metadata$condition)))
stopifnot(!anyNA(parsed))
design <- model.matrix(~ donor + condition, parsed)
stopifnot(ncol(design) == 4L, qr(design)$rank == 4L)

# Two numerically equivalent donor labels must still form distinct blocks.
write.csv(metadata[1:4, ], sheet, row.names = FALSE, quote = FALSE)
paired <- read_sample_metadata(sheet, metadata$sample_id[1:4])
stopifnot(nlevels(factor(paired$donor)) == 2L)
stopifnot(ncol(model.matrix(~ donor + condition, paired)) == 3L)

failure <- tryCatch(read_sample_metadata(sheet, c("missing", "1", "002", "2")), error = identity)
stopifnot(inherits(failure, "error"))
unlink(sheet)
message("Passed: literal sample/donor/condition labels, pairing, order and count alignment.")
