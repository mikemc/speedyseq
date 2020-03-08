# Attribution: The documentation and coding strategy for `tax_glom()` and
# `tip_glom()` are modified from phyloseq
# https://github.com/joey711/phyloseq/blob/master/R/transform_filter-methods.

#' Agglomerate taxa of the same type.
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
#' This is the speedyseq reimplementation of phyloseq's tax_glom function. It
#' is designed to produce identical results but this has not been thoroughly
#' tested. Please report any discrepancies.
#'
#' Documentation and general strategy derived from `phyloseq::tax_glom()`.
#'
#' @usage tax_glom(physeq, taxrank=rank_names(physeq)[1], NArm=TRUE,
#'   bad_empty=c(NA, "", " ", "\t"))
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} or
#'   \code{\link{tax_table}}. (NOTE: Currently only works on phyloseq objects)
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
#' 
#' @return A taxonomically-agglomerated, optionally-pruned, object with class
#' matching the class of \code{physeq}.
#'
#' @seealso
#' \code{\link{phyloseq::tax_glom}}
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
#' # data(GlobalPatterns)
#' # ## print the available taxonomic ranks
#' # colnames(tax_table(GlobalPatterns))
#' # ## agglomerate at the Family taxonomic rank
#' # (x1 <- tax_glom(GlobalPatterns, taxrank="Family") )
#' # ## How many taxa before/after agglomeration?
#' # ntaxa(GlobalPatterns); ntaxa(x1)
#' # ## Look at enterotype dataset...
#' # data(enterotype)
#' # ## print the available taxonomic ranks. Shows only 1 rank available, not useful for tax_glom
#' # colnames(tax_table(enterotype))
tax_glom <- function(physeq, 
                     taxrank = rank_names(physeq)[1],
                     NArm = TRUE, 
                     bad_empty = c(NA, "", " ", "\t")) {
  if (is.null(access(physeq, "tax_table"))) {
    stop("The tax_glom() function requires that physeq contain a taxonomyTable")
  }
  if (!taxrank[1] %in% rank_names(physeq)) {
    stop("Bad taxrank argument. Must be among the values of rank_names(physeq)")
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
  physeq <- merge_taxa_vec(physeq, tax_strings, tax_adjust = 0L)
  # "Empty" the taxonomy values to the right of the rank, using NA_character_.
  if (rank_idx < length(rank_names(physeq))) {
    bad_ranks <- seq(rank_idx + 1, length(rank_names(physeq)))
    tax_table(physeq)[, bad_ranks] <- NA_character_
  }
  physeq
}

tax_glom0 <- function(physeq, 
                     taxrank = rank_names(physeq)[1],
                     NArm = TRUE, 
                     bad_empty = c(NA, "", " ", "\t")) {
  #### This part is identical to phyloseq's tax_glom
	# Error if tax_table slot is empty
	if( is.null(access(physeq, "tax_table")) ){
		stop("The tax_glom() function requires that physeq contain a taxonomyTable")
	}
	# Error if bad taxrank
	if( !taxrank[1] %in% rank_names(physeq) ){
		stop("Bad taxrank argument. Must be among the values of rank_names(physeq)")
	}
	# Make a vector from the taxonomic data.
	CN  <- which( rank_names(physeq) %in% taxrank[1] )
	tax <- as(access(physeq, "tax_table"), "matrix")[, CN]
	# if NArm is TRUE, remove the empty, white-space, NA values from 
	if( NArm ){
		keep_species <- names(tax)[ !(tax %in% bad_empty) ]
		physeq <- prune_taxa(keep_species, physeq)
	}
	# Concatenate data up to the taxrank column, use this for agglomeration
	tax <- as(access(physeq, "tax_table"), "matrix")[, 1:CN, drop=FALSE]
	tax <- apply(tax, 1, function(i){paste(i, sep=";_;", collapse=";_;")})
  #### **Speedyseq changes start here**
  # Note, `tax` is a named vector whose names are the OTU/taxa names and whose
  # elements are strings of the ranks up to `taxrank` such as
  # "Archaea;_;Crenarchaeota;_;Thaumarchaeota;_;Cenarchaeales" 
  ## Make the new OTU table
  # Work with taxa as rows
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
    # Note that we need to flip back to samples as rows at the end
    needs_flip <- TRUE
  } else {
    needs_flip <- FALSE
  }
  otu <- otu_table(physeq)
  # Compute new table with base::rowsum()
  new_otu <- otu_table(rowsum(otu, tax), taxa_are_rows = TRUE)
  # The archetype is the most abundant OTU in each group. The first OTU is
  # chosen in the case of ties, which should be consistent with original
  # phyloseq behavior.
  archetypes <- vapply(
    split(taxa_sums(otu), tax),
    function (x) names(x)[which.max(x)],
    "a"
  )
  stopifnot(all.equal(names(archetypes), rownames(new_otu)))
  taxa_names(new_otu) <- unname(archetypes)
  ## Make the new phyloseq object
  # Replacing the original otu_table with the new, smaller table will
  # automatically prune the taxonomy, tree, and refseq to the smaller set of
  # archetypal OTUs
  otu_table(physeq) <- new_otu
  # "Empty" the taxonomy values to the right of the rank, using NA_character_.
	if (CN < length(rank_names(physeq))) {
    bad_ranks <- seq(CN + 1, length(rank_names(physeq)))
    tax_table(physeq)[, bad_ranks] <- NA_character_
  }
	## Return.
  if (needs_flip) {
    physeq <- t(physeq)
  }
	return(physeq)
}

#' Agglomerate closely-related taxa using phylogeny-derived distances
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
#' Documentation and general strategy derived from `phyloseq::tip_glom()`.
#'
#' Speedyseq note: \code{\link[stats]{hclust}} is faster than the default
#' `hcfun`; set `method = "average"` to get equivalent clustering.
#'
#' @param physeq (Required). A \code{\link{phyloseq-class}}, containing a
#'   phylogenetic tree. Alternatively, a phylogenetic tree
#'   \code{\link[ape]{phylo}} will also work. (NOTE: currently only works on
#'   phyloseq objects)
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
  if (! "phyloseq" %in% class(physeq))
    stop('Currently only implemented for phyloseq objects. Use `phyloseq::tip_glom()` for phylo objects.')
  tree <- access(physeq, "phy_tree")
  if (is.null(tree)) {
    stop("The tip_glom() function requires that `physeq` contain a phylogenetic tree")
  }
  d <- castor::get_all_pairwise_distances(
    tree,
    only_clades = taxa_names(tree),
    check_input = FALSE
  )
  rownames(d) <- colnames(d) <- taxa_names(tree)
  d <- as.dist(d)
  psclust <- cutree(as.hclust(hcfun(d, ...)), h = h)
  merge_taxa_vec(physeq, psclust, tax_adjust = tax_adjust)
}
