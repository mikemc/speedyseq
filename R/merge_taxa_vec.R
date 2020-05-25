# merge_taxa_vec --------------------------------------------------------------

#' Merge taxa in groups (vectorized version)
#' 
#' Merge taxa in `x` into a smaller set of taxa defined by the vector `group`.
#' Taxa whose value in `group` is NA will be dropped. New taxa will be named
#' according to the most abundant taxon in each group (`phyloseq` and
#' `otu_table` objects) or the first taxon in each group (all other phyloseq
#' component objects).
#'
#' The `tax_adjust` argument controls the handling of taxonomic disagreements
#' within groups. Setting `tax_adjust == 0` causes no adjustment; the taxonomy
#' of the new group is set to the archetype taxon (see below). Otherwise,
#' disagreements within a group at a given rank cause the values at lower ranks
#' to be set to `NA`. If `tax_adjust == 1` (the default), then a rank where all
#' taxa in the group are already NA is not counted as a disagreement, and lower
#' ranks may be kept if the taxa agree. This corresponds to the original
#' phyloseq behavior. If `tax_adjust == 2`, then these NAs are treated as a
#' disagreement; all ranks are set to NA after the first disagreement or NA.
#' 
#' @param x A phyloseq object or component object 
#' @param group a vector with one element for each taxon in `physeq` that
#'   defines the new groups. see `base::rowsum()`.
#' @param tax_adjust 0: no adjustment; 1: phyloseq-compatible adjustment; 2:
#'   conservative adjustment (see Details)
#' @export
#'
#' @seealso
#' \code{\link[phyloseq]{merge_taxa}}
#'
#' \code{\link{tax_glom}}
#'
#' @rdname merge_taxa_vec-methods
#' @examples
#' \dontrun{
#' # Create 97% OTUs from a phyloseq object with reference sequences
#' # Assume ps is a phyloseq object with refseq slot
#' dna <- refseq(ps)
#' nproc <- 1 # Increase to use multiple processors
#' aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
#' d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
#' clusters <- DECIPHER::IdClusters(
#'   d, 
#'   method = "complete",
#'   cutoff = 0.03,
#'   processors = nproc
#' )
#' ps0 <- merge_taxa_vec(
#'   ps, 
#'   group = clusters$cluster,
#'   tax_adjust = 2
#' )
#' }
setGeneric("merge_taxa_vec", 
  function(x, 
           group, 
           tax_adjust = 1L)
    standardGeneric("merge_taxa_vec")
)

#' @rdname merge_taxa_vec-methods
setMethod("merge_taxa_vec", "phyloseq",
  function(x, group, tax_adjust = 1L) {
    stopifnot(ntaxa(x) == length(group))
    stopifnot(tax_adjust %in% c(0L, 1L, 2L))
    # drop taxa with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding taxa will be dropped")
      x <- prune_taxa(!is.na(group), x)
      group <- group[!is.na(group)]
    }
    # Get the merged otu table with new taxa named by most abundant
    otu <- merge_taxa_vec(otu_table(x), group)
    # Adjust taxonomy if necessary
    if (!is.null(x@tax_table) & tax_adjust != 0) {
      tax <- merge_taxa_vec(tax_table(x), group, tax_adjust = tax_adjust)
      # Taxa in `tax` are in same order as in `otu` but are named by first in
      # group and so need to be renamed
      taxa_names(tax) <- taxa_names(otu)
    }
    # Create the new phyloseq object. Replacing the original otu_table with
    # the new, smaller table will automatically prune the taxonomy, tree, and
    # refseq to the smaller set of archetypal taxa.
    otu_table(x) <- otu
    if (exists("tax"))
      tax_table(x) <- tax
    x
  }
)

#' @rdname merge_taxa_vec-methods
setMethod("merge_taxa_vec", "otu_table",
  function(x, group) {
    stopifnot(ntaxa(x) == length(group))
    # Work with taxa as rows
    if (!taxa_are_rows(x)) {
      x <- t(x)
      needs_flip <- TRUE
    } else {
      needs_flip <- FALSE
    }
    # drop taxa with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding taxa will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    # New taxa names are the most abundant taxon in each group; in the case of
    # ties, the first taxon is chosen. Original group order is maintained.
    new_names <- data.table::data.table(
      taxon = taxa_names(x), 
      sum = taxa_sums(x),
      group = group
    ) %>%
      .[, by = group, .(archetype = taxon[which.max(sum)])]
    # Compute new table with base::rowsum(). The call to rowsum() makes the
    # rownames the group names; reorder = FALSE keeps the taxa in the same
    # order as `new_names`.
    otu <- otu_table(rowsum(x, group, reorder = FALSE), taxa_are_rows = TRUE)
    if (needs_flip) {
      otu <- t(otu)
    }
    stopifnot(all.equal(as.character(new_names$group), taxa_names(otu)))
    taxa_names(otu) <- new_names$archetype
    otu
  }
)

#' @rdname merge_taxa_vec-methods
setMethod("merge_taxa_vec", "taxonomyTable",
  function(x, group, tax_adjust = 1L) {
    stopifnot(ntaxa(x) == length(group))
    # drop taxa with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding taxa will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    if (tax_adjust == 0L)
      return(merge_taxa_vec_pseudo(x, group))
    else if (tax_adjust == 1L)
      na_bad <- FALSE
    else if (tax_adjust == 2L)
      na_bad <- TRUE
    k <- length(rank_names(x))
    # bad_string is used to temporarily mark bad values in the tax table
    bad_string <- paste0("BAD", Sys.time())
    tax <- x %>% 
      as("matrix") %>%
      data.table::as.data.table(keep.rownames = "taxon") %>%
      .[, group := group] %>%
      .[,
        # Reduce to one row per group; compute new names and bad ranks
        by = group, .SDcols = rank_names(x),
        c(
          # New taxa names are the first taxon in each group
          .(taxon = taxon[1]), 
          lapply(.SD, bad_or_unique, bad = bad_string)
        )
      ] %>%
      .[, !"group"] %>%
      # Propagate bad ranks downwards and convert to NAs
      tibble::column_to_rownames("taxon") %>%
      apply(1, bad_flush_right, bad = bad_string, na_bad = na_bad, k = k) %>% 
      t %>%
      tax_table
    return(tax)
  }
)

#' @rdname merge_taxa_vec-methods
setMethod("merge_taxa_vec", "phylo", 
  function(x, group) {merge_taxa_vec_pseudo(x, group)}
)

#' @rdname merge_taxa_vec-methods
setMethod("merge_taxa_vec", "XStringSet",
  function(x, group) {merge_taxa_vec_pseudo(x, group)}
)

#' Pseudo taxa merging for phylo and XStringSet objects
#'
#' Returns `x` pruned to the first taxon of each group defined in `group`.
#'
#' @param x a phylo or XStringSet object
#' @param group a vector with one element for each taxon in `physeq` that
#'   defines the new groups
#' @keywords internal
merge_taxa_vec_pseudo <- function(x, group) {
  stopifnot(ntaxa(x) == length(group))
  # drop taxa with `is.na(group)`
  if (anyNA(group)) {
    warning("`group` has missing values; corresponding taxa will be dropped")
    x <- prune_taxa(!is.na(group), x)
    group <- group[!is.na(group)]
  }
  # Archetypes are the first taxon in each group
  archetypes <- data.table::data.table(taxon = taxa_names(x), group = group) %>%
    .[, by = group, .(taxon = taxon[1])] %>%
    .$taxon
  prune_taxa(archetypes, x)
}

# helper functions ------------------------------------------------------------

#' Reduce a vector x to its unique value or the value of `bad`
#'
#' Helper for `merge_taxa_vec()`
#'
#' @param x a vector
#' @param bad the string representing a bad value
#' @keywords internal
bad_or_unique <- function(x, bad = "BAD") {
  if (length(unique(x)) == 1)
    x[[1]]
  else
    bad
}

#' Replace all values with NA upon seeing a bad value
#'
#' Helper for `merge_taxa_vec()`
#'
#' @param x a vector
#' @param bad the string representing a bad value
#' @param na_bad whether NAs should also be treated as bad
#' @param k the index to which values are flushed
#' @keywords internal
bad_flush_right <- function(x, bad = "BAD", na_bad = FALSE, k = length(x)) {
  if (na_bad) 
    which_bad <- which(x == bad | is.na(x))
  else 
    which_bad <- which(x == bad)
  if (length(which_bad) > 0)
    x[seq(min(which_bad), k)] <- NA
  x
}
