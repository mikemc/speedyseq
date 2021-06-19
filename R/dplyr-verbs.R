# filter ----------------------------------------------------------------------

#' Subset taxa using values in the taxonomy table
#'
#' This function makes it possible to choose a subset of taxa using values in
#' the taxonomy table as in `dplyr::filter()`. The taxonomy table is converted
#' to a tibble using `ps_tibble()` and the arguments in `...` are passed
#' directly to `dplyr::filter()`. The taxonomy table in the phyloseq object is
#' then updated to contain just the subset taxa.
#' 
#' @param x A `phyloseq` or `taxonomyTable` object
#' @param ... Expressions passed to `dplyr::filter()`
#' 
#' @export
setGeneric("filter_tax_table", 
  function(x, ...) standardGeneric("filter_tax_table")
)

#' @rdname filter_tax_table
setMethod("filter_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% filter_tax_table(...)
    x
  })

#' @rdname filter_tax_table
setMethod("filter_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::filter(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' Subset samples using values in the sample data
#'
#' This function is a wrapper around `dplyr::filter()` that provides a
#' convenient way to subset samples using the `sample_data(x)`.
#' 
#' @param x A `phyloseq` or `sample_data` object
#' @param ... Expressions passed to `dplyr::filter()`
#' 
#' @export
setGeneric("filter_sample_data", 
  function(x, ...) standardGeneric("filter_sample_data")
)

#' @rdname filter_sample_data
setMethod("filter_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% filter_sample_data(...)
    x
  })

#' @rdname filter_sample_data
setMethod("filter_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::filter(...) %>%
      {suppressMessages(sample_data(.))}
  })

# mutate ----------------------------------------------------------------------

#' Create, modify, and delete columns in the taxonomy table
#'
#' This function is a wrapper around `dplyr::mutate()` that provides a
#' convenient way to modify `tax_table(x)`.
#' 
#' @param x A `phyloseq` or `taxonomyTable` object
#' @param ... Expressions passed to `dplyr::mutate()`
#' 
#' @export
setGeneric("mutate_tax_table", 
  function(x, ...) standardGeneric("mutate_tax_table")
)

#' @rdname mutate_tax_table
setMethod("mutate_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% mutate_tax_table(...)
    x
  })

#' @rdname mutate_tax_table
setMethod("mutate_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::mutate(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' Create, modify, and delete columns in the sample data
#'
#' This function is a wrapper around `dplyr::mutate()` that provides a
#' convenient way to modify `sample_data(x)`.
#' 
#' @param x A `phyloseq` or `sample_data` object
#' @param ... Expressions passed to `dplyr::mutate()`
#' 
#' @export
setGeneric("mutate_sample_data", 
  function(x, ...) standardGeneric("mutate_sample_data")
)

#' @rdname mutate_sample_data
setMethod("mutate_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% mutate_sample_data(...)
    x
  })

#' @rdname mutate_sample_data
setMethod("mutate_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::mutate(...) %>%
      {suppressMessages(sample_data(.))}
  })
