#' @importFrom tibble as_tibble
#' @export
tibble::as_tibble

#' Coerce phyloseq objects to tibbles (tbl_df objects)
#'
#' These functions extend the `as_tibble()` function defined by the tibble
#' package to work on phyloseq objects or phyloseq-component objects.
#'
#' @param x a phyloseq object or component
#' @param pivot whether to pivot the data frame to long format
#' @param tax whether to include taxonomy data
#'
#' @name as_tibble-speedyseq
NULL

#' @export
#' @describeIn as_tibble-speedyseq method for otu_table objects
as_tibble.otu_table <- function(x, pivot = TRUE, .name_repair = base::make.unique) {
  mat <- x %>%
    as("matrix")
  if (pivot) {
    nms <- c(
      getOption("speedyseq.tibble_sample"), 
      getOption("speedyseq.tibble_otu"), 
      getOption("speedyseq.tibble_abundance")
    )
    if (anyDuplicated(nms)) 
      stop("OTU, sample, and abundance column names must be unique")
    if (!taxa_are_rows(x))
      mat <- t(mat)
    tb <- mat %>%
      tibble::as_tibble(rownames = getOption("speedyseq.tibble_otu")) %>%
      tidyr::pivot_longer(
        cols = sample_names(x),
        names_to = getOption("speedyseq.tibble_sample"), 
        values_to = getOption("speedyseq.tibble_abundance"), 
      )
  } else {
    rn <- ifelse(taxa_are_rows(x), 
      getOption("speedyseq.tibble_otu"), 
      getOption("speedyseq.tibble_sample")
    )
    tb <- mat %>%
      tibble::as_tibble(rownames = rn) %>%
      rlang::set_names(., 
        vctrs::vec_as_names(names(.), repair = .name_repair)
      )
  }
  tb
}

#' @export
#' @describeIn as_tibble-speedyseq method for sample_data objects
as_tibble.sample_data <- function(x, .name_repair = base::make.unique) {
  x %>%
    as("data.frame") %>%
    tibble::as_tibble(rownames = getOption("speedyseq.tibble_sample")) %>%
    rlang::set_names(., 
      vctrs::vec_as_names(names(.), repair = .name_repair)
    )
}

#' @export
#' @describeIn as_tibble-speedyseq method for taxonomyTable objects
as_tibble.taxonomyTable <- function(x, .name_repair = base::make.unique) {
  x %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = getOption("speedyseq.tibble_otu")) %>%
    rlang::set_names(., 
      vctrs::vec_as_names(names(.), repair = .name_repair)
    )
}

#' @export
#' @describeIn as_tibble-speedyseq method for DNAStringSet objects
as_tibble.XStringSet <- function(x) {
  nms <- c(
    getOption("speedyseq.tibble_otu"), 
    getOption("speedyseq.tibble_sequence")
  )
  if (anyDuplicated(nms)) 
    stop("OTU and sequence column names must be unique")
  x %>%
    as.character %>%
    tibble::enframe(nms[1], nms[2])
}

#' @export
#' @describeIn as_tibble-speedyseq method for phyloseq objects
as_tibble.phyloseq <- function(x, 
                               tax = TRUE, 
                               ref = FALSE, 
                               .name_repair = base::make.unique) {
  nms <- c(
    getOption("speedyseq.tibble_otu"), 
    getOption("speedyseq.tibble_sample"), 
    getOption("speedyseq.tibble_abundance")
  )
  if (anyDuplicated(nms)) 
    stop("OTU, sample, and abundance column names must be unique")
  tb <- otu_table(x) %>% as_tibble(pivot = TRUE)
  # Add sample data if it exists
  sam <- access(x, "sam_data")
  if (!is.null(sam)) {
    tb <- tb %>% 
      dplyr::left_join(
        sam %>% as_tibble(.name_repair = .name_repair),
        by = nms[2],
        suffix = c("", ".sam")
      )
  }
  # Add tax_table if it exists and tax = TRUE
  tt <- access(x, "tax_table")
  if (tax & !is.null(tt)) {
    tb <- tb %>% 
      dplyr::left_join(
        tt %>% as_tibble(.name_repair = .name_repair),
        by = nms[1], 
        suffix = c("", ".tax")
      )
  }
  # Add refseq if it exists and ref = TRUE
  rs <- access(x, "refseq")
  if (ref & !is.null(rs)) {
    tb <- tb %>% 
      dplyr::left_join(
        rs %>% as_tibble,
        by = nms[1], 
        suffix = c("", ".ref")
      )
  }
  tb
}

# #' @export
#> setGeneric("as_tibble")

#' @keywords internal
as_tibble <- function(x, ...) {
  UseMethod("as_tibble")
}

setMethod("as_tibble", "otu_table", as_tibble.otu_table)
setMethod("as_tibble", "sample_data", as_tibble.sample_data)
setMethod("as_tibble", "taxonomyTable", as_tibble.taxonomyTable)
setMethod("as_tibble", "XStringSet", as_tibble.XStringSet)
setMethod("as_tibble", "phyloseq", as_tibble.phyloseq)
