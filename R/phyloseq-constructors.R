#' Construct phyloseq objects from tibbles
#'
#' These methods extend phyloseq's constructor functions to construct phyloseq
#' components from tibbles (objects with class "tbl_df"). 
#'
#' Since tibbles cannot have row names, the sample or taxon identifiers must be
#' contained in a regular column. Speedyseq currently always uses the first
#' column for the identifiers that would normally be taken from the row names
#' by phyloseq's built-in constructors. Thus the first column is assumed to
#' contain the sample names for `sample_data()` and the OTU/taxa names for
#' `tax_table()`; for `otu_table()`, the first column is assumed to contain the
#' sample names if `taxa_are_rows = TRUE` and the taxa names if `taxa_are_rows
#' = FALSE`, with the other identifier being taken from the remaining column
#' names.
#' 
#' @param object A tibble whose first column contains the sample or taxa ids
#' @param taxa_are_rows Logical; `TRUE` if rows correspond to taxa and `FALSE`
#' if rows correspond to samples
#'
#' @name tibble-constructors
#' 
#' @seealso 
#' [`tibble::tbl_df`]
#' [phyloseq::otu_table()]
#' [phyloseq::sample_data()]
#' [phyloseq::tax_table()]
#'
#' @examples
#' \dontrun{
#' # Read a .csv file with readr, which creates an object of class `tbl_df`
#' tbl <- readr::read_csv("path/to/otu_table.csv")
#' # Inspect and check if taxa are rows and that the first column contains the
#' # sample names or the taxa/OTU names
#' head(tbl) 
#' # Create a phyloseq `otu_table` object
#' otu <- otu_table(tbl, taxa_are_rows = FALSE)
#'
#' # Read a .csv file with readr, which creates an object of class `tbl_df`
#' tbl <- readr::read_csv("path/to/sample_data.csv")
#' # Inspect and check that the first column contains the sample names
#' head(tbl) 
#' # Create a phyloseq `sample_data` object
#' sam <- sample_data(tbl)
#' }

#' @rdname tibble-constructors
#' @export
setMethod("otu_table", "tbl_df", function (object, taxa_are_rows) {
  if (taxa_are_rows) {
    message(paste0(
      "Assuming first column, `", names(object)[1], 
      "`, contains the taxa names"))
  } else {
    message(paste0(
      "Assuming first column, `", names(object)[1], 
      "`, contains the sample names"))
  }
  object %>%
    tibble::column_to_rownames(var = names(object)[1]) %>%
    as("matrix") %>%
    otu_table(taxa_are_rows)
})

#' @rdname tibble-constructors
#' @export
setMethod("sample_data", "tbl_df", function(object) {
  message(paste0(
    "Assuming first column, `", names(object)[1], 
    "`, contains the sample names"))
  object %>%
    tibble::column_to_rownames(var = names(object)[1]) %>%
    sample_data
})

#' @rdname tibble-constructors
#' @export
setMethod("tax_table", "tbl_df", function (object) {
  message(paste0(
    "Assuming first column, `", names(object)[1], 
    "`, contains the taxa names"))
  object %>%
    tibble::column_to_rownames(var = names(object)[1]) %>%
    as("matrix") %>%
    tax_table
})

