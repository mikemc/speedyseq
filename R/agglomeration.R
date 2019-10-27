# tax glom --------------------------------------------------------------------

#' Agglomerate taxa of the same type.
#'
#' This method merges species that have the same taxonomy at a certain
#' taxonomic rank.  Its approach is analogous to \code{\link{tip_glom}}, but
#' uses categorical data instead of a tree. In principal, other categorical
#' data known for all taxa could also be used in place of taxonomy, but for the
#' moment, this must be stored in the \code{taxonomyTable} of the data. Also,
#' columns/ranks to the right of the rank chosen to use for agglomeration will
#' be replaced with \code{NA}, because they should be meaningless following
#' agglomeration.
#'
#' NOTE: This is the speedyseq reimplementation of phyloseq's tax_glom
#' function. It is designed to produce identical results but this has not been
#' thoroughly tested. Please report any discrepancies!
#'
#' @usage tax_glom(physeq, taxrank=rank_names(physeq)[1], NArm=TRUE,
#' bad_empty=c(NA, "", " ", "\t"))
#'
#' @param physeq (Required). \code{\link{phyloseq-class}} or \code{\link{otu_table}}.
#'
#' @param taxrank A character string specifying the taxonomic level
#'  that you want to agglomerate over.
#'  Should be among the results of \code{rank_names(physeq)}.
#'  The default value is \code{rank_names(physeq)[1]},
#'  which may agglomerate too broadly for a given experiment.
#'  You are strongly encouraged to try different values for this argument.
#'
#' @param NArm (Optional). Logical, length equal to one. Default is \code{TRUE}.
#'  CAUTION. The decision to prune (or not) taxa for which you lack categorical
#'  data could have a large effect on downstream analysis. You may want to
#'  re-compute your analysis under both conditions, or at least think carefully
#'  about what the effect might be and the reasons explaining the absence of 
#'  information for certain taxa. In the case of taxonomy, it is often a result 
#'  of imprecision in taxonomic designation based on short phylogenetic sequences
#'  and a patchy system of nomenclature. If this seems to be an issue for your
#'  analysis, think about also trying the nomenclature-agnostic \code{\link{tip_glom}}
#'  method if you have a phylogenetic tree available.
#'
#' @param bad_empty (Optional). Character vector. Default: \code{c(NA, "", " ", "\t")}.
#'  Defines the bad/empty values 
#'  that should be ignored and/or considered unknown. They will be removed
#'  from the internal agglomeration vector derived from the argument to \code{tax},
#'  and therefore agglomeration will not combine taxa according to the presence
#'  of these values in \code{tax}. Furthermore, the corresponding taxa can be
#'  optionally pruned from the output if \code{NArm} is set to \code{TRUE}.
#' 
#' @return A taxonomically-agglomerated, optionally-pruned, object with class matching
#' the class of \code{physeq}.
#'
#' @seealso
#' \code{\link{tip_glom}}
#' 
#' \code{\link{prune_taxa}}
#' 
#' \code{\link{merge_taxa}}
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
tax_glom <- function(physeq, taxrank=rank_names(physeq)[1],
					NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){
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
    ## Make the new OTU table
    # Work with taxa as rows
    if (!taxa_are_rows(physeq)) {
        physeq <- phyloseq::t(physeq)
        # Note that we need to flip back to samples as rows at the end
        needs_flip <- TRUE
    } else {
        needs_flip <- FALSE
    }
    # Starting point is a tibble with rows as taxa, to be able to combine taxa
    # with the dplyr::summarize_*() functions
    otu <- otu_table(physeq)
    tb <- otu %>%
        as("matrix") %>%
        tibble::as_tibble(rownames = "OTU")
    # We want to name each new taxon (group of merged OTUs) by its "archetype",
    # the most abundant OTU in the group
    tb <- tb %>%
        tibble::add_column(Tax = tax, Sum = taxa_sums(physeq)) %>%
        dplyr::group_by(Tax)
    # Name new taxa by the most abundant OTU; pick the first OTU in case of
    # ties (to be consistent with phyloseq)
    new_taxa_names <- tb %>% 
        dplyr::top_n(1, Sum) %>%
        dplyr::slice(1) %>%
        dplyr::select(Tax, OTU)
    # Sum abundances and rename taxa
    tb0 <- tb %>%
        dplyr::summarize_at(dplyr::vars(sample_names(physeq)), sum) %>%
        dplyr::left_join(new_taxa_names, by = "Tax") %>%
        dplyr::select(OTU, dplyr::everything(), -Tax)
    # Put back into phyloseq form
    mat <- tb0 %>%
        dplyr::select(-OTU) %>%
        as("matrix")
    rownames(mat) <- tb0$OTU
    otu0 <- otu_table(mat, taxa_are_rows = TRUE)
    ## Make the new phyloseq object
    # Replacing the original otu_table with the new, smaller table will
    # automatically prune the taxonomy, tree, and refseq to the smaller set of
    # archetypal otus
    otu_table(physeq) <- otu0
	# "Empty" the taxonomy values to the right of the rank, using
    # NA_character_.
	if (CN < length(rank_names(physeq))) {
        bad_ranks <- seq(CN + 1, length(rank_names(physeq)))
        tax_table(physeq)[, bad_ranks] <- NA_character_
    }
    if (needs_flip) {
        physeq <- phyloseq::t(physeq)
    }
	physeq
}

# tree glom --------------------------------------------------------------------

#' Agglomerate taxa to a given phylogenetic resolution.
#'
#' This method merges taxa that are sufficiently closely related according to
#' `phy_tree(physeq)`. It uses
#' \code{\link{castor::collapse_tree_at_resolution}} to determine the groups of
#' taxa that are to be aggregated and to produce a new tree with each group
#' collapsed to a single tip. The new taxa are named according to the most
#' abundant taxon of each group. The reference sequences and taxonomy table
#' reflect these "archetypal" taxa.
#'
#' This method is similar to `phyloseq::tip_glom()` except that it uses the
#' tree directly, rather than as the source of a distance matrix that is used
#' for hierarchical clustering. The clustering step of `tip_glom()` is slow and
#' memory intensive when `ntaxa(physeq)` is large. Since clustering does not
#' need to be performed, `tree_glom()` is much more efficient.
#'
#' @param physeq \code{\link{phyloseq-class}}.
#' @param resolution Phylogenetic resolution at which to collapse the tree
#' @param shorten Whether to shorten the terminal branches to the collapsed
#'   nodes
#' @param criterion Criterion to determine whether a node is collapsed
#'
#' @return A taxonomically-agglomerated, optionally-pruned, object with class matching
#' the class of \code{physeq}.
#'
#' @seealso
#' \code{\link{castor::collapse_tree_at_resolution}} for more information about
#' tree-collapsing method and parameters 
#'
#' \code{\link{phyloseq::tip_glom}}
#'
#' \code{\link{tax_glom}}
#' 
#' @export
#'
#' @examples
#' data(GlobalPatterns)
#' ntaxa(GlobalPatterns)
#' ps <- tree_glom(GlobalPatterns, 0.03)
#' ntaxa(ps)
#'
tree_glom <- function(physeq, resolution, shorten = TRUE,
    criterion = "max_tip_depth") {

    # Find the groups of taxa to collapse -------------------------------------

    tree <- phy_tree(physeq)
    collapsed <- castor::collapse_tree_at_resolution(tree, resolution, 
        by_edge_count = FALSE, shorten = shorten,
        rename_collapsed_nodes = TRUE, criterion = criterion)

    # collapsed$collapsed_nodes contains the nodes that were collapsed --- the
    # tips from each of these nodes tells us the groups of taxa that we need to
    # merge
    collapsed_groups <- collapsed$collapsed_nodes %>%
        purrr::map(~castor::get_subtree_at_node(tree, .)$subtree$tip.label)
    # Name each group by the member that labels the tip in the collapsed tree
    tips <- collapsed_groups %>%
        purrr::map_chr(~ .[. %in% collapsed$tree$tip.label])
    names(collapsed_groups) <- tips
    tips <- tibble::enframe(collapsed_groups, "collapsed_tip", "taxon") %>%
        tidyr::unnest("taxon")
    # This is for the collapsed nodes only; for the other taxa, we simply have
    # collapsed_tip = taxon
    tips <- taxa_names(physeq) %>%
        setdiff(tips$taxon) %>%
        rlang::set_names() %>%
        tibble::enframe("collapsed_tip", "taxon") %>%
        dplyr::bind_rows(tips)
    # Now we have the groups, and the name in the collapsed tree

    # Archetypes --------------------------------------------------------------

    # We want to name each new taxon (group of merged OTUs) by its "archetype",
    # the most abundant OTU in the group

    archetypes <- taxa_sums(physeq) %>%
        tibble::enframe("taxon", "sum") %>%
        dplyr::left_join(tips, by = "taxon") %>%
        dplyr::group_by(collapsed_tip) %>%
        dplyr::mutate(rank = rank(-sum, ties.method = "first")) %>%
        dplyr::filter(rank == 1) %>%
        dplyr::select(collapsed_tip, archetype = taxon)

    tips <- dplyr::left_join(tips, archetypes, by = "collapsed_tip")

    # New tree ----------------------------------------------------------------

    # If shorten == TRUE, then we should keep the new collapsed tree, and
    # rename the tips to be the archetypes
    if (shorten == TRUE) {
        new_tree <- collapsed$tree
        new_tree$tip.label <- dplyr::left_join(
                tibble::tibble(collapsed_tip = collapsed$tree$tip.label),
                archetypes, by = "collapsed_tip") %>%
            dplyr::pull(archetype)
    }

    # If shorten == FALSE, then we should stick to the original tree, but prune
    # to the archetypes
    if (shorten == FALSE) 
        new_tree <- prune_taxa(archetypes$archetype, tree)

    # New OTU table -----------------------------------------------------------

    # Starting point is a tibble with rows as taxa, to be able to combine taxa
    # with the dplyr::summarize_*() functions
    if (!taxa_are_rows(physeq)) {
        physeq <- phyloseq::t(physeq)
        # Note that we need to flip back to samples as rows at the end
        needs_flip <- TRUE
    } else {
        needs_flip <- FALSE
    }
    otutb <- otu_table(physeq) %>%
        as("matrix") %>%
        tibble::as_tibble(rownames = "taxon") %>%
        dplyr::left_join(tips, by = "taxon") %>%
        dplyr::group_by(archetype) %>%
        dplyr::summarize_at(sample_names(physeq), sum)
    # Put back into phyloseq form
    mat <- otutb %>%
        dplyr::select(-archetype) %>%
        as("matrix")
    rownames(mat) <- otutb$archetype
    new_otu <- otu_table(mat, taxa_are_rows = TRUE)

    # New tax table -----------------------------------------------------------

    # TODO: make a new tax table manually, to figure out which tax ranks are
    # still valid for the new OTU.



    # Phyloseq object ---------------------------------------------------------

    # Replacing the original otu_table with the new, smaller table will
    # automatically prune the taxonomy, tree, and refseq to the smaller set of
    # archetypal otus
    phy_tree(physeq) <- new_tree
    otu_table(physeq) <- new_otu
    if (needs_flip) {
        physeq <- phyloseq::t(physeq)
    }
	physeq
}

