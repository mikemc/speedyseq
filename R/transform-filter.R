# Function titles copied from phyloseq's documentation

# Transform sample counts -----------------------------------------------------

#' Transform abundance data in an `otu_table`, sample-by-sample
#' 
#' Version of [phyloseq::transform_sample_counts()] that allows a
#' purrr-style anonymous function.
#'
#' This function simply calls [purrr::as_mapper()] on `.f` and passes
#' the resulting function on to
#' [phyloseq::transform_sample_counts()].
#'
#' @param physeq [phyloseq-class()] or [ape::phylo()].
#' @param .f A function or formula that can be converted to a function by
#'   [purrr::as_mapper()]
#' @param ... Additional arguments passed on to `.f`
#' 
#' @seealso 
#' [phyloseq::transform_sample_counts()]
#' [purrr::as_mapper()]
#'
#' @export
#' @examples
#' data(GlobalPatterns)
#' # Filter low prevalence taxa, then convert to proportions
#' gp.prop <- GlobalPatterns %>%
#'   filter_taxa2(~ sum(. > 0) > 5) %>%
#'   transform_sample_counts(~ . / sum(.))
transform_sample_counts <- function(physeq, .f, ...) {
  fun <- purrr::as_mapper(.f)
  phyloseq::transform_sample_counts(physeq, fun, ...)
}

# Filter taxa -----------------------------------------------------------------

#' Filter taxa based on across-sample OTU abundance criteria
#' 
#' Variations of [phyloseq::filter_taxa()] that allows a purrr-style
#' anonymous function.
#'
#' `filter_taxa()` simply calls [purrr::as_mapper()] on `.f` and
#' passes the resulting function on to [phyloseq::filter_taxa()].
#' `filter_taxa2()` also sets `prune = TRUE`, which is convenient when passing
#' a phyloseq object in a pipe chain (see example).
#'
#' @param physeq [phyloseq-class()] or [ape::phylo()].
#' @param .f A function or formula that can be converted to a function by
#'   [purrr::as_mapper()]
#' @param prune A logical. If `FALSE`, then this function returns a logical
#'   vector specifying the taxa that passed the filter; if `TRUE`, then this
#'   function returns the pruned phyloseq object.
#'
#' @export
#' @seealso 
#' [phyloseq::filter_taxa()]
#' [purrr::as_mapper()]
#' 
#' @examples
#' data(GlobalPatterns)
#' # Filter low prevalence taxa and then convert to proportions
#' gp.prop <- GlobalPatterns %>%
#'   filter_taxa2(~ sum(. > 0) > 5) %>%
#'   transform_sample_counts(~ . / sum(.))
filter_taxa <- function(physeq, .f, prune = FALSE) {
  fun <- purrr::as_mapper(.f)
  phyloseq::filter_taxa(physeq, fun, prune = prune)
}

#' @export
#' @describeIn filter_taxa Sets `prune = TRUE`
filter_taxa2 <- function(physeq, .f) {
  fun <- purrr::as_mapper(.f)
  phyloseq::filter_taxa(physeq, fun, prune = TRUE)
}

