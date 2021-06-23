speedyseq_default_options <- list(
  speedyseq.psmelt_class = "data.frame"
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(speedyseq_default_options) %in% names(op))
  if (any(toset)) 
    options(speedyseq_default_options[toset])
  invisible()
}
