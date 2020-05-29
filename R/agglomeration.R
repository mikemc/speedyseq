# Attribution: The documentation and coding strategy for `tax_glom()` and
# `tip_glom()` are modified from phyloseq
# https://github.com/joey711/phyloseq/blob/master/R/transform_filter-methods.

#' Agglomerate taxa using taxonomy.
#'
#' This method merges species that have the same taxonomy at a certain
#' taxonomic rank. Its approach is analogous to \code{\link{tip_glom}}, but
#' uses categorical data instead of a tree. In principal, other categorical
#' data known for all taxa could also be used in place of taxonomy, but for the
#' moment, this must be stored in the \code{taxonomyTable} of the data. Also,
#' columns/ranks to the right of the rank chosen to use for agglomeration will
#' be replaced with \code{NA}, because they should be meaningless following
#' agglomeration.
#'
#' This is the speedyseq reimplementation of `phyloseq::tax_glom()`. It should
#' produce results that are identical to phyloseq up to taxon order. 
#'
#' If `x` is a phyloseq object with a phylogenetic tree, then the new taxa will
#' be ordered as they are in the tree. Otherwise, the taxa order can be
#' controlled by the `reorder` argument, which behaves like the `reorder`
#' argument in \code{\link[base]{rowsum}}. `reorder = FALSE` will keep taxa in
#' the original order determined by when the member of each group first appears
#' in `taxa_names(x)`; `reorder = TRUE` will order new taxa alphabetically
#' according to taxonomy (string of concatenated rank values).
#'
#' Acknowledgements: Documentation and general strategy derived from
#' `phyloseq::tax_glom()`.
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} or
#'   \code{\link{tax_table}}.
#' @param taxrank A character string specifying the taxonomic level that you
#'   want to agglomerate over. Should be among the results of
#'   \code{rank_names(physeq)}. The default value is
#'   \code{rank_names(physeq)[1]}, which may agglomerate too broadly for a
#'   given experiment. You are strongly encouraged to try different values for
#'   this argument.
#' @param NArm (Optional). Logical, length equal to one. Default is
#'   \code{TRUE}.  CAUTION. The decision to prune (or not) taxa for which you
#'   lack categorical data could have a large effect on downstream analysis.
#'   You may want to re-compute your analysis under both conditions, or at
#'   least think carefully about what the effect might be and the reasons
#'   explaining the absence of information for certain taxa. In the case of
#'   taxonomy, it is often a result of imprecision in taxonomic designation
#'   based on short phylogenetic sequences and a patchy system of nomenclature.
#'   If this seems to be an issue for your analysis, think about also trying
#'   the nomenclature-agnostic \code{\link{tip_glom}} method if you have a
#'   phylogenetic tree available.
#' @param bad_empty (Optional). Character vector. Default: \code{c(NA, "", " ",
#'   "\t")}. Defines the bad/empty values that should be ignored and/or
#'   considered unknown. They will be removed from the internal agglomeration
#'   vector derived from the argument to \code{tax}, and therefore
#'   agglomeration will not combine taxa according to the presence of these
#'   values in \code{tax}. Furthermore, the corresponding taxa can be
#'   optionally pruned from the output if \code{NArm} is set to \code{TRUE}.
#' @param reorder Logical specifying whether to reorder the taxa by taxonomy
#'   strings or keep initial order. Ignored if `physeq` has a phylogenetic
#'   tree.
#' 
#' @return A taxonomically-agglomerated, optionally-pruned, object with class
#' matching the class of \code{physeq}.
#'
#' @seealso
#' \code{\link[phyloseq]{tax_glom}}
#' 
#' \code{\link{tip_glom}}
#' 
#' \code{\link{prune_taxa}}
#' 
#' \code{\link{merge_taxa_vec}}
#' 
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' # print the available taxonomic ranks
#' colnames(tax_table(GlobalPatterns))
#' # agglomerate at the Family taxonomic rank
#' (x1 <- tax_glom(GlobalPatterns, taxrank="Family"))
#' # How many taxa before/after agglomeration?
#' ntaxa(GlobalPatterns); ntaxa(x1)
tax_glom <- function(physeq, 
                     taxrank = rank_names(physeq)[1],
                     NArm = TRUE, 
                     bad_empty = c(NA, "", " ", "\t"),
                     reorder = FALSE) {
  if (is.null(access(physeq, "tax_table"))) {
    stop("`tax_glom()` requires that `physeq` contain a taxonomy table")
  }
  if (!taxrank[1] %in% rank_names(physeq)) {
    stop("Bad `taxrank` argument. Must be among the values of `rank_names(physeq)`")
  }
  rank_idx <- which(rank_names(physeq) %in% taxrank[1])
  # if NArm is TRUE, remove taxa whose value for taxrank is in bad_empty
  if (NArm) {
    bad_taxa <- tax_table(physeq)[, rank_idx] %in% bad_empty
    physeq <- prune_taxa(!bad_taxa, physeq)
  }
  # For each taxon, make a string containing its full taxonomy, which will
  # define the groups to be merged
  tax_strings <- apply(
    tax_table(physeq)[, 1:rank_idx, drop = FALSE],
    1, 
    function(x) {paste(x, collapse=";")}
  )
  # Merge taxa with speedyseq's vectorized merging function
  physeq <- merge_taxa_vec(physeq, tax_strings, tax_adjust = 0L, 
    reorder = reorder)
  # "Empty" the taxonomy values to the right of the rank, using NA_character_.
  if (rank_idx < length(rank_names(physeq))) {
    bad_ranks <- seq(rank_idx + 1, length(rank_names(physeq)))
    tax_table(physeq)[, bad_ranks] <- NA_character_
  }
  physeq
}

#' Agglomerate taxa using phylogeny-derived distances.
#' 
#' All tips of the tree separated by a phylogenetic cophenetic distance smaller
#' than \code{h} will be agglomerated into one taxon.
#' 
#' Can be used to create a non-trivial OTU Table, if a phylogenetic tree is
#' available.
#'
#' By default, simple ``greedy'' single-linkage clustering is used. It is
#' possible to specify different clustering approaches by setting \code{hcfun}
#' and its parameters in `...`. In particular, complete-linkage clustering
#' appears to be used more commonly for OTU clustering applications.
#'
#' The merged taxon is named according to the "archetype" defined as the the
#' most abundant taxon (having the largest value of `taxa_sums(physeq)`. The
#' tree and refseq objects are pruned to the archetype taxa.
#'
#' Speedyseq note: \code{\link[stats]{hclust}} is faster than the default
#' `hcfun`; set `method = "average"` to get equivalent clustering.
#'
#' Acknowledgements: Documentation and general strategy derived from
#' `phyloseq::tip_glom()`.
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}}, containing a
#'   phylogenetic tree. Alternatively, a phylogenetic tree
#'   \code{\link[ape]{phylo}} will also work.
#' @param h (Optional). Numeric scalar of the height where the tree should be
#'   cut. This refers to the tree resulting from hierarchical clustering of the
#'   distance matrix, not the original phylogenetic tree. Default value is
#'   \code{0.2}.
#' @param hcfun (Optional). A function. The (agglomerative, hierarchical)
#'   clustering function to use. The default is \code{\link[cluster]{agnes}}
#'   for phyloseq compatiblity.
#' @param tax_adjust 0: no adjustment; 1: phyloseq-compatible adjustment; 2:
#'   conservative adjustment (see \code{\link{merge_taxa_vec}} for details)
#' @param ... (Optional). Additional named arguments to pass to \code{hcfun}. 
#' @return An instance of the \code{\link{phyloseq-class}}.  Or alternatively,
#'   a \code{\link{phylo}} object if the \code{physeq} argument was just a
#'   tree.  In the expected-use case, the number of OTUs will be fewer (see
#'   \code{\link{ntaxa}}), after merging OTUs that are related enough to be
#'   called the same OTU. 
#'
#' @seealso 
#'
#' \code{\link[phyloseq]{tip_glom}}
#'
#' \code{\link{tree_glom}} for direct phylogenetic merging
#' 
#' \code{\link{merge_taxa_vec}}
#' 
#' \code{\link[cluster]{agnes}}
#' 
#' \code{\link[stats]{hclust}}
#' 
#' \code{\link[castor]{get_all_pairwise_distances}}
#' 
#' \code{\link[ape]{phylo}}
#'
#' @export
#'
#' @examples 
#' data("esophagus")
#' esophagus <- prune_taxa(taxa_names(esophagus)[1:25], esophagus)
#' plot_tree(esophagus, label.tips="taxa_names", size="abundance", 
#'   title="Before tip_glom()")
#' plot_tree(tip_glom(esophagus, h=0.2), label.tips="taxa_names", 
#'   size="abundance", title="After tip_glom()")
#' 
#' # *speedyseq only:* Demonstration of different `tax_adjust` behaviors
#' data(GlobalPatterns)
#' 
#' set.seed(20190421)
#' ps <- prune_taxa(sample(taxa_names(GlobalPatterns), 2e2), GlobalPatterns)
#'
#' ps1 <- tip_glom(ps, 0.1, tax_adjust = 1)
#' ps2 <- tip_glom(ps, 0.1, tax_adjust = 2)
#' tax_table(ps1)[c(108, 136, 45),]
#' tax_table(ps2)[c(108, 136, 45),]
tip_glom <- function(physeq, 
                     h = 0.2, 
                     # hcfun = hclust, # Need method = "average" for compat
                     hcfun = cluster::agnes,
                     tax_adjust = 1L,
                     ...) {
  tree <- access(physeq, "phy_tree")
  if (is.null(tree)) {
    stop("`tip_glom()` requires that physeq contain a phylogenetic tree")
  }
  d <- castor::get_all_pairwise_distances(
    tree,
    only_clades = taxa_names(tree),
    check_input = FALSE
  )
  rownames(d) <- colnames(d) <- taxa_names(tree)
  d <- stats::as.dist(d)
  psclust <- stats::cutree(stats::as.hclust(hcfun(d, ...)), h = h)
  merge_taxa_vec(physeq, psclust, tax_adjust = tax_adjust)
}

#' Agglomerate taxa at a given phylogenetic resolution.
#'
#' This function merges taxa that are sufficiently closely related in the tree
#' in `phy_tree(physeq)`. It is similar to `tip_glom()`, except that it uses
#' the tree directly, rather than as the basis for hierarchical clustering.
#'
#' This function uses \code{\link[castor]{collapse_tree_at_resolution}} from
#' the castor package to determine the groups of taxa to be collapsed; see this
#' function for details about the `resolution` and `criterion` parameters.
#'
#' New taxa are named according to the most abundant taxon of each group. The
#' tree, reference sequences, and (by default) the taxonomy table reflect these
#' "archetype" taxa.
#'
#' @param physeq \code{\link{phyloseq-class}} or \code{\link[ape]{phylo}}.
#' @param resolution Phylogenetic resolution at which to merge taxa. 
#' @param criterion Criterion for determining whether to collapse taxa. See
#'   \code{\link[castor]{collapse_tree_at_resolution}} for details.
#' @param tax_adjust 0: no adjustment; 1: phyloseq-compatible adjustment; 2:
#'   conservative adjustment (see \code{\link{merge_taxa_vec}} for details)
#'
#' @return A \code{physeq} object with new taxa reflecting the
#'   phylogenetically merged groups.
#'
#' @seealso
#' \code{\link[castor]{collapse_tree_at_resolution}} for more information about
#'   the `resolution` and `criterion` parameters.
#'
#' \code{\link{merge_taxa_vec}} for more about `tax_adjust` and general merging
#' behavior
#'
#' \code{\link{tip_glom}} for indirect phylogenetic merging
#'
#' \code{\link{tax_glom}}
#' 
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' ps1 <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
#' ntaxa(ps1)
#' ps2 <- tree_glom(ps1, 0.05)
#' ntaxa(ps2)
#'
#' library(dplyr)
#' library(ggtree)
#' library(cowplot)
#'
#' plot1 <- phy_tree(ps1) %>%
#'   ggtree +
#'   geom_tiplab() +
#'   geom_label(aes(x = branch, label = round(branch.length, 4)))
#' plot2 <- phy_tree(ps2) %>%
#'   ggtree +
#'   geom_tiplab() +
#'   geom_label(aes(x = branch, label = round(branch.length, 4)))
#' plot_grid(plot1, plot2)
tree_glom <- function(physeq, 
                      resolution, 
                      criterion = "max_tip_depth",
                      tax_adjust = 1L) {
  tree <- access(physeq, "phy_tree")
  if (is.null(tree)) {
    stop("`tree_glom()` requires that `physeq` contain a phylogenetic tree")
  }

  # Find the groups of taxa to collapse ---------------------------------------
  # Collapse the tree (results include info about the collapsed groups)
  collapsed <- castor::collapse_tree_at_resolution(
    tree, 
    resolution, 
    by_edge_count = FALSE, 
    shorten = FALSE,
    # `rename_collapsed_nodes = TRUE` can give errors
    rename_collapsed_nodes = FALSE,
    criterion = criterion
  )

  # If no collapsed nodes, return physeq as is
  if (length(collapsed$collapsed_nodes) == 0)
    return(physeq)

  # Data frame mapping collapsed tips to group names
  collapsed_tbl <- tibble::tibble(
    collapsed_node = collapsed$collapsed_nodes
    ) %>%
    dplyr::mutate(
      collapsed_tip = purrr::map(collapsed_node,
        ~castor::get_subtree_at_node(tree, .)$subtree$tip.label
      )
    ) %>%
    tidyr::unnest(collapsed_tip) %>%
    # Name the group by the first taxon, rather than the node name, to ensure no
    # name conflicts. Archetypes are computed in merge_taxa_vec later on.
    dplyr::group_by(collapsed_node) %>%
    dplyr::mutate(
      group = utils::head(collapsed_tip, 1)
    ) %>%
    dplyr::ungroup()

  # Compute merged ps object --------------------------------------------------
  # Will use merge_taxa_vec with a `group` vector reflecting the
  # phylogenetically similar groups; this will automatically prune the
  # reference sequences and tree to the archetypes
  group <- tibble::tibble(taxon = taxa_names(physeq)) %>%
    dplyr::left_join(collapsed_tbl, 
      by = c("taxon" = "collapsed_tip")) %>%
    dplyr::mutate(group = ifelse(is.na(group), taxon, group)) %>%
    dplyr::select(taxon, group) %>%
    tibble::deframe()
  merge_taxa_vec(physeq, group, tax_adjust = tax_adjust)
}
