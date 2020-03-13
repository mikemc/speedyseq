# merge_taxa_vec --------------------------------------------------------------

#' Merge taxa in groups (vectorized version)
#' 
#' Merge taxa in `x` into a smaller set of taxa defined by the vector `group`.
#' Taxa whose value in `group` is NA will be dropped. New taxa will be named
#' either by the most abundant taxon in each group (`name == "max"`).
#'
#' The `tax_adjust` argument controls the handling of taxonomic disagreements
#' within groups. Setting `tax_adjust == 0` causes no adjustment; the taxonomy
#' of the new group is set to the archetype taxon (see below). Setting
#' `tax_adjust == 1` corresponds to the original phyloseq behavior.
#' Disagreements within a group at any rank cause the values at lower ranks to
#' be set to `NA`, but ranks where all taxa in the group are already NA are not
#' treated as disagreements. As a result, taxonomic assignments where e.g. the
#' family is NA but the genus is not are be propagated to the new table.
#' Setting `tax_adjust == 2` (the default) goes further, setting all ranks to
#' NA after the first disagreement or NA.
#' 
#' @param x A phyloseq object or component object 
#' @param group A vector with one element for each taxon in `physeq` that
#'   defines the new groups. See `base::rowsum()`.
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
    stopifnot(is.vector(group))
    stopifnot(ntaxa(x) == length(group))
    stopifnot(tax_adjust %in% c(0L, 1L, 2L))
    # Get the merged otu table with new taxa set to the archetype (max) names
    new_otu_table <- merge_taxa_vec(otu_table(x), group)
    # Create a new taxonomy if requested. Otherwise the archetypes' taxonomies
    # will be used.
    if (!is.null(x@tax_table) & tax_adjust != 0L) {
      if (tax_adjust == 1L)
        na_bad <- FALSE
      else if (tax_adjust == 2L)
        na_bad <- TRUE
      k <- length(rank_names(x))
      # need a string not in the tax table to mark bad values
      bad_string <- paste0("BAD", Sys.time())
      new_tax_mat <- tax_table(x)@.Data %>% 
        as.data.frame(stringsAsFactors = FALSE) %>%
        stats::aggregate(
          by = list(group = group),
          bad_or_unique, bad = bad_string
        ) %>%
        tibble::column_to_rownames("group") %>%
        apply(1, bad_flush_right, bad = bad_string, na_bad = na_bad, k = k) %>%
        t
      rownames(new_tax_mat) <- taxa_names(new_otu_table)
    }
    # Create the new phyloseq object. Replacing the original otu_table with
    # the new, smaller table will automatically prune the taxonomy, tree, and
    # refseq to the smaller set of archetypal OTUs
    otu_table(x) <- new_otu_table
    if (exists("new_tax_mat"))
      tax_table(x) <- tax_table(new_tax_mat)
    x
  }
)

#' @rdname merge_taxa_vec-methods
setMethod("merge_taxa_vec", "otu_table",
  function(x, group) {
    stopifnot(is.vector(group))
    stopifnot(ntaxa(x) == length(group))
    # Work with taxa as rows
    if (!taxa_are_rows(x)) {
      otu <- t(x)
      needs_flip <- TRUE
    } else {
      otu <- x
      needs_flip <- FALSE
    }
    # drop taxa with `is.na(group)`
    if (any(is.na(group))) {
      warning("`group` contains NAs; corresponding taxa will be dropped")
      otu <- otu[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    # Get list of new taxa names from the abundant taxon in each group. In
    # the case of ties, the first taxon is chosen
    new_taxa_names <- split(taxa_sums(otu), group) %>%
      vapply(function(y) names(y)[which.max(y)], "a")
    # Compute new table with base::rowsum(). The call to rowsum() makes the
    # rownames the group names; reorder = TRUE puts the taxa in the same order
    # as `new_taxa_names`.
    otu <- otu_table(rowsum(otu, group, reorder = TRUE), taxa_are_rows = TRUE)
    if (needs_flip) {
      otu <- t(otu)
    }
    stopifnot(all.equal(names(new_taxa_names), taxa_names(otu)))
    taxa_names(otu) <- unname(new_taxa_names)
    otu
  }
)

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
