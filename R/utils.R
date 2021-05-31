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

# orient_taxa -----------------------------------------------------------------

#' Orient a phyloseq object or otu table to have taxa as rows or as columns
#'
#' Puts the phyloseq or otu-table object `x` in the orientation (taxa as rows
#' or as columns) specified by `as`. This is useful when passing the otu table
#' on to functions that require the abundance matrix to have a specific
#' orientation and are unaware of the `taxa_are_rows(x)` property.
#'
#' @param x A phyloseq or otu-table object
#' @param as The matrix dimension that is desired for taxa. Must be
#'   "rows" for rows and "columns" or "cols" for columns.
#'
#' @export
#'
#' @examples
#' data(soilrep)
#' taxa_are_rows(soilrep)
#' x <- soilrep %>% orient_taxa(as = "columns")
#' taxa_are_rows(x)
setGeneric("orient_taxa", 
  function(x, as) standardGeneric("orient_taxa")
)

orient_taxa_default <- function(x, as) {
  stopifnot(as %in% c("rows", "columns", "cols"))
  if (identical(as, "rows")) {
    if (!taxa_are_rows(x))
      x <- t(x)
  } else {
    if (taxa_are_rows(x))
      x <- t(x)
  }
  x
}

#' @rdname orient_taxa
setMethod("orient_taxa", "otu_table", orient_taxa_default)

#' @rdname orient_taxa
setMethod("orient_taxa", "phyloseq", orient_taxa_default)
