# Attribution: The documentation for `psmelt()` and the first 54 lines of the
# function's code are from phyloseq,
# https://github.com/joey711/phyloseq/blob/master/R/plot-methods.R

#' Melt phyloseq data object into large data.frame
#'
#' The psmelt function is a specialized melt function for melting phyloseq
#' objects (instances of the phyloseq class), usually for producing graphics
#' with [ggplot2::ggplot2]. The naming conventions used in downstream phyloseq
#' graphics functions have reserved the following variable names that should
#' not be used as the names of [sample_variables()] or taxonomic
#' [rank_names()]. These reserved names are `c("Sample", "Abundance", "OTU")`.
#' Also, you should not have identical names for sample variables and taxonomic
#' ranks. That is, the intersection of the output of the following two
#' functions [sample_variables()], [rank_names()] should be an empty vector
#' (e.g.  `intersect(sample_variables(physeq), rank_names(physeq))`). All of
#' these potential name collisions are checked-for and renamed automatically
#' with a warning. However, if you (re)name your variables accordingly ahead of
#' time, it will reduce confusion and eliminate the warnings.
#'
#' The `as` argument allows choosing between three classes for tabular data:
#'
#' - [data.frame-class] is chosen by "data.frame" or "df"
#' - [data.table-class] is chosen by "data.table" or "dt"
#' - [tbl_df-class] is chosen by "tbl_df", "tbl", or "tibble"
#'
#' The default is "data.frame" and can be overridden by setting the
#' "speedyseq.psmelt_class" option.
#'
#' Note that "melted" phyloseq data is stored much less efficiently, and so RAM
#' storage issues could arise with a smaller dataset (smaller number of
#' samples/OTUs/variables) than one might otherwise expect.  For common sizes
#' of graphics-ready datasets, however, this should not be a problem.  Because
#' the number of OTU entries has a large effect on the RAM requirement, methods
#' to reduce the number of separate OTU entries -- for instance by
#' agglomerating OTUs based on phylogenetic distance using [tip_glom()] -- can
#' help alleviate RAM usage problems. This function is made user-accessible for
#' flexibility, but is also used extensively by plot functions in phyloseq.
#'
#' @param physeq An [otu_table-class] or [phyloseq-class]; most useful for
#'   phyloseq-class.
#' @param as Class of the output table; see Details.
#'
#' @return A table of the specified class
#'
#' @seealso
#' [plot_bar()]
#' 
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
#' mdf = psmelt(gp.ch)
#' nrow(mdf)
#' ncol(mdf)
#' colnames(mdf)
#' head(rownames(mdf))
#' # Create a ggplot similar to
#' library("ggplot2")
#' p = ggplot(mdf, aes(x=SampleType, y=Abundance, fill=Genus))
#' p = p + geom_bar(color="black", stat="identity", position="stack")
#' print(p)
psmelt <- function(physeq, as = getOption("speedyseq.psmelt_class")) {
  stopifnot(as %in% c("data.frame", "df", "data.table", "dt", "tbl_df", "tbl",
      "tibble"))
  # Access covariate names from object, if present
  if(!inherits(physeq, "phyloseq")){
    rankNames = NULL
    sampleVars = NULL
  } else {
    # Still might be NULL, but attempt access
    rankNames = rank_names(physeq, FALSE)
    sampleVars = sample_variables(physeq, FALSE) 
  }
  # Define reserved names
  reservedVarnames = c("Sample", "Abundance", "OTU")  
  # type-1a conflict: between sample_data 
  # and reserved psmelt variable names
  type1aconflict = intersect(reservedVarnames, sampleVars)
  if(length(type1aconflict) > 0){
    wh1a = which(sampleVars %in% type1aconflict)
    new1a = paste0("sample_", sampleVars[wh1a])
    # First warn about the change
    warning("The sample variables: \n",
            paste(sampleVars[wh1a], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new1a, collapse=", "), "\n",
            "to avoid conflicts with special phyloseq plot attribute names.")
    # Rename the sample variables.
    colnames(sample_data(physeq))[wh1a] <- new1a
  }
  # type-1b conflict: between tax_table
  # and reserved psmelt variable names
  type1bconflict = intersect(reservedVarnames, rankNames)
  if(length(type1bconflict) > 0){
    wh1b = which(rankNames %in% type1bconflict)
    new1b = paste0("taxa_", rankNames[wh1b])
    # First warn about the change
    warning("The rank names: \n",
            paste(rankNames[wh1b], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new1b, collapse=", "), "\n",
            "to avoid conflicts with special phyloseq plot attribute names.")
    # Rename the conflicting taxonomic ranks
    colnames(tax_table(physeq))[wh1b] <- new1b
  }
  # type-2 conflict: internal between tax_table and sample_data
  type2conflict = intersect(sampleVars, rankNames)
  if(length(type2conflict) > 0){
    wh2 = which(sampleVars %in% type2conflict)
    new2 = paste0("sample_", sampleVars[wh2])
    # First warn about the change
    warning("The sample variables: \n",
            paste0(sampleVars[wh2], collapse=", "), 
            "\n have been renamed to: \n",
            paste0(new2, collapse=", "), "\n",
            "to avoid conflicts with taxonomic rank names.")
    # Rename the sample variables
    colnames(sample_data(physeq))[wh2] <- new2
  }
  # Enforce OTU table orientation. Redundant-looking step
  # supports "naked" otu_table as `physeq` input.
  otutab = otu_table(physeq)
  if(!taxa_are_rows(otutab)){otutab <- t(otutab)}
  ## Speedyseq changes start here
  # Convert the otu table to a tibble in tall form (one sample-taxon obsevation
  # per row)
  dtb <- otutab %>% 
    as("matrix") %>%
    data.table::as.data.table(keep.rownames = "OTU") %>%
    data.table::melt(id.vars = c("OTU"), variable.name = "Sample", 
      value.name = "Abundance", variable.factor = FALSE)
  # Add the sample data if it exists
  if (!is.null(sampleVars)) {
    sam <- sample_data(physeq) %>%
      as("data.frame") %>% 
      data.table::as.data.table(keep.rownames = "Sample")
    dtb <- dtb[sam, on = .(Sample = Sample)]
  }
  # Add the tax table if it exists
  if (!is.null(rankNames)) {
    tax <- tax_table(physeq) %>%
      as("matrix") %>%
      data.table::as.data.table(keep.rownames = "OTU")
    dtb <- dtb[tax, on = .(OTU = OTU)]
  }
  # Arrange by Abundance, then OTU names (to approx. phyloseq behavior)
  dtb <- dtb %>%
    data.table::setorder(-Abundance, OTU)
  # Return as requested class
  if (as %in% c("data.table", "dt")) {
    dtb
  } else if (as %in% c("data.frame", "df")) {
    dtb %>% as("data.frame")
  } else if (as %in% c("tbl_df", "tbl", "tibble")) {
    dtb <- dtb %>% tibble::as_tibble()
    attr(dtb, ".internal.selfref") <- NULL
    dtb
  } else {
    stop("Invalid output class specified by `as`.")
  }
}
