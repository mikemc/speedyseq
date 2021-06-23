# vim: foldmethod=marker

# Some of the documentation for these functions is modified from the
# corresponding dplyr functions (https://github.com/tidyverse/dplyr), MIT
# license RStudio and others.

#' Mutating joins of the taxonomy table or sample data
#'
#' @description
#' `r lifecycle::badge("experimental")`
#'
#' These functions are wrappers around dplyr joining operations (see
#' `dplyr::mutate-joins`) whose first (`x`) argument is either the taxonomy
#' table or sample data.
#'
#' @details
#' When joining by taxonomy table, the .otu column name can be used to match
#' taxa names.
#' When joining by sample data, the .sample column name can be used to match
#' sample names.
#' 
#' @param x A `phyloseq`, `taxonomyTable`, or `sample_data` object
#' @param ... Arguments passed to respective dplyr joining operation
#'
#' @return An object of the same type as `x`, with added columns
#'
#' @seealso 
#' [dplyr::left_join()]
#' [dplyr::inner_join()]
#'
#' @name join-phyloseq
#'
#' @examples
#' data(GlobalPatterns)
#' 
#' GlobalPatterns %>% sample_variables
#' ps1 <- GlobalPatterns %>%
#'   select_sample_data(!contains("Barcode"))
#' y <- GlobalPatterns %>%
#'   sample_data %>%
#'   select_sample_data(contains("Barcode")) %>%
#'   ps_tibble
#' ps2 <- ps1 %>% left_join_sample_data(y, by = ".sample")
#' ps2 %>% sample_variables
NULL

# left_join ----------------------------------------------------------------{{{

#' @rdname join-phyloseq
#' @export
setGeneric("left_join_tax_table", 
  function(x, ...) standardGeneric("left_join_tax_table")
)

setMethod("left_join_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% left_join_tax_table(...)
    x
  })

setMethod("left_join_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::left_join(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname join-phyloseq
#' @export
setGeneric("left_join_sample_data", 
  function(x, ...) standardGeneric("left_join_sample_data")
)

setMethod("left_join_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% left_join_sample_data(...)
    x
  })

setMethod("left_join_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::left_join(...) %>%
      {suppressMessages(sample_data(.))}
  })

# }}}

# inner_join ---------------------------------------------------------------{{{

#' @rdname join-phyloseq
#' @export
setGeneric("inner_join_tax_table", 
  function(x, ...) standardGeneric("inner_join_tax_table")
)

setMethod("inner_join_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% inner_join_tax_table(...)
    x
  })

setMethod("inner_join_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::inner_join(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname join-phyloseq
#' @export
setGeneric("inner_join_sample_data", 
  function(x, ...) standardGeneric("inner_join_sample_data")
)

setMethod("inner_join_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% inner_join_sample_data(...)
    x
  })

setMethod("inner_join_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::inner_join(...) %>%
      {suppressMessages(sample_data(.))}
  })

# }}}
