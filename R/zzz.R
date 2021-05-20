speedyseq_default_options <- list(
  speedyseq.psmelt_class = "data.frame",
  # column names for as_tibble methods; match psmelt by default
  speedyseq.tibble_sample = ".sample",
  speedyseq.tibble_otu = ".otu",
  speedyseq.tibble_abundance = ".abundance",
  speedyseq.tibble_sequence = ".sequence"
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(speedyseq_default_options) %in% names(op))
  if (any(toset)) 
    options(speedyseq_default_options[toset])
  invisible()
}
