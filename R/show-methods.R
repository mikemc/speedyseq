# show/print ------------------------------------------------------------------

# Overwrite show method for phyloseq objects

#' @export
setMethod("show", "phyloseq", function(object){
	cat("phyloseq-class experiment-level object", fill = TRUE)
	# print otu_table (always there).
	cat(paste("otu_table()   OTU Table:         [ ", ntaxa(otu_table(object)), " taxa and ", nsamples(otu_table(object)), " samples ]", sep = ""), fill = TRUE)	
	# print Sample Data if there
	if(!is.null(sample_data(object, FALSE))){
        cat(paste("sample_data() Sample Data:       [ ", dim(sample_data(object))[1], " samples by ", 
	        dim(sample_data(object))[2], 
            " sample variables ]", sep = ""), fill = TRUE)
	}
	# print tax table if there	
	if(!is.null(tax_table(object, FALSE))){
        cat(paste("tax_table()   Taxonomy Table:    [ ", dim(tax_table(object))[1], " taxa by ", 
	        dim(tax_table(object))[2], 
            " taxonomic ranks ]", sep = ""), fill = TRUE)
	}
	# print tree if there
	if(!is.null(phy_tree(object, FALSE))){
        cat(paste("phy_tree()    Phylogenetic Tree: [ ", ntaxa(phy_tree(object)), " tips and ", 
	        phy_tree(object)$Nnode,
            " internal nodes ]", sep = ""),
        	fill = TRUE
        ) 
	}
	# print refseq summary if there
	if(!is.null(refseq(object, FALSE))){
        cat(paste("refseq()      ", class(refseq(object))[1], ":      [ ", ntaxa(refseq(object)), " reference sequences ]", sep = ""), fill=TRUE)
	}
  if (taxa_are_rows(object))
    cat("taxa are rows", fill = TRUE)
  else
    cat("taxa are columns", fill = TRUE)
})


# glimpse ---------------------------------------------------------------------

# Fix output of tibble::glimpse() on sample_data and taxonomyTabl eobjects by
# first converting to a data.frame. Note, currently the sample and taxa names
# will not be shown.

#' @export
glimpse.sample_data <- function(x, width = NULL, ...) {
  x %>% 
    as("data.frame") %>% 
    tibble::glimpse(width = width, ...)
}

#' @export
glimpse.taxonomyTable <- function(x, width = NULL, ...) {
  x %>% 
    as.data.frame %>% 
    tibble::glimpse(width = width, ...)
}
