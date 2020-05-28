# select_taxa -----------------------------------------------------------------

#' Select a subset of taxa in a specified order where possible
#'
#' Select (a subset of) taxa; if `x` allows taxa to be reordered, then taxa are
#' given in the specified order.
#'
#' This is a simple selector function that is like `prune_taxa(taxa, x)` when
#' `taxa` is a character vector but always gives the taxa in the order `taxa`
#' if possible (that is, except for phy_tree's and phyloseq objects that
#' contain phy_tree's).
#'
#' @param x A phyloseq object or phyloseq component object
#' @param taxa Character vector of taxa to select, in requested order
#' @param reorder Logical specifying whether to use the order in `taxa` (TRUE)
#'   or keep the order in `taxa_names(x)` (FALSE)
#' @keywords internal
#' @rdname select_taxa-methods
setGeneric("select_taxa", 
  function(x, taxa, reorder = TRUE) standardGeneric("select_taxa")
)

#' @rdname select_taxa-methods
setMethod("select_taxa", signature("sample_data", "character"), 
  function(x, taxa) {
    stopifnot(!anyDuplicated(taxa))
    x
  }
)

#' @rdname select_taxa-methods
setMethod("select_taxa", signature("otu_table", "character"), 
  function(x, taxa, reorder = TRUE){
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder)
      taxa <- intersect(taxa_names(x), taxa)
    if (taxa_are_rows(x)) {
      x[taxa, , drop=FALSE]
    } else {
      x[, taxa, drop=FALSE]
    }
  }
)

#' @rdname select_taxa-methods
setMethod("select_taxa", signature("taxonomyTable", "character"), 
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder)
      taxa <- intersect(taxa_names(x), taxa)
		x[taxa, , drop=FALSE]
	}
)

#' @rdname select_taxa-methods
setMethod("select_taxa", signature("XStringSet", "character"), 
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder)
      taxa <- intersect(taxa_names(x), taxa)
		x[taxa]
	}
)

#' @rdname select_taxa-methods
setMethod("select_taxa", signature("phylo", "character"), 
  function(x, taxa) {
    # NOTE: `reorder` argument silently ignored if supplied
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    ape::keep.tip(x, taxa)
  }
)

#' @rdname select_taxa-methods
setMethod("select_taxa", signature("phyloseq", "character"), 
  function(x, taxa, reorder = TRUE) {
    stopifnot(!anyDuplicated(taxa))
    stopifnot(all(taxa %in% taxa_names(x)))
    if (!reorder)
      taxa <- intersect(taxa_names(x), taxa)
    otu_table(x) <- select_taxa(otu_table(x), taxa)
    phyloseq:::index_reorder(x, index_type = "taxa")
  }
)
