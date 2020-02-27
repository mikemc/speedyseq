# Attribution: The documentation for `tax_glom()` is from phyloseq,
# https://github.com/joey711/phyloseq/blob/master/R/transform_filter-methods.

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
tax_glom <- function(physeq, 
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
