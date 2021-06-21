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

#' Create, modify, and delete columns in the taxonomy table or sample data
#'
#' These functions are wrappers around `dplyr::mutate()` and
#' `dplyr::transmute()` that provide convenient ways to modify `tax_table(x)`
#' and `sample_data(x)` as well as the sample and taxa names.
#'
#' When modifying the taxonomy table, the `.otu` column name can be used to set
#' new taxa names.
#' When modifying the sample data, the `.sample` column name can be used to set
#' new taxa names.
#'
#' The experimental arguments to `dplyr::mutate()` of `.keep`, `.before`, and
#' `.after` are not currently supported and may result in errors or unexpected
#' behavior.
#' 
#' @param x A `phyloseq`, `taxonomyTable`, or `sample_data` object
#' @param ... Expressions passed to `dplyr::mutate()` or `dplyr::transmute()`
#'
#' @name mutate-phyloseq
#'
#' @examples
#' data(GlobalPatterns)
#' 
#' ps <- GlobalPatterns %>%
#'   transmute_sample_data(SampleType, sample_sum = sample_sums(.)) %>%
#'   filter_sample_data(SampleType %in% c("Feces", "Skin", "Tongue")) %>%
#'   filter_tax_table(Kingdom == "Bacteria") %>%
#'   tax_glom("Phylum") %>%
#'   transmute_tax_table(Kingdom, Phylum, .otu = Phylum)
#' sample_data(ps)
#' tax_table(ps)
NULL

## mutate tax table

#' @rdname mutate-phyloseq
#' @export
setGeneric("mutate_tax_table", 
  function(x, ...) standardGeneric("mutate_tax_table")
)

#' @rdname mutate-phyloseq
setMethod("mutate_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::mutate(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname mutate-phyloseq
setMethod("mutate_tax_table", "phyloseq",
  function(x, ...) {
    # Mutating the tax table may create new taxa names, which need to get
    # updated in the phyloseq object prior to updating its tax table
    new_tax <- tax_table(x) %>% mutate_tax_table(...)
    taxa_names(x) <- taxa_names(new_tax)
    tax_table(x) <- new_tax
    x
  })

## transmute tax table

#' @rdname mutate-phyloseq
#' @export
setGeneric("transmute_tax_table", 
  function(x, ...) standardGeneric("transmute_tax_table")
)

#' @rdname mutate-phyloseq
setMethod("transmute_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      # Make sure the .otu column is always included; new names can still be
      # set if .otu occurs in the ... expressions
      dplyr::transmute(.otu, ...) %>%
      dplyr::relocate(.otu) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname mutate-phyloseq
setMethod("transmute_tax_table", "phyloseq",
  function(x, ...) {
    # Mutating the tax table may create new taxa names, which need to get
    # updated in the phyloseq object prior to updating its tax table
    new_tax <- tax_table(x) %>% transmute_tax_table(...)
    taxa_names(x) <- taxa_names(new_tax)
    tax_table(x) <- new_tax
    x
  })

## mutate sample data

#' @rdname mutate-phyloseq
#' @export
setGeneric("mutate_sample_data", 
  function(x, ...) standardGeneric("mutate_sample_data")
)

#' @rdname mutate-phyloseq
setMethod("mutate_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::mutate(...) %>%
      {suppressMessages(sample_data(.))}
  })

#' @rdname mutate-phyloseq
setMethod("mutate_sample_data", "phyloseq",
  function(x, ...) {
    # Mutating the sample data may create new sample names, which need to get
    # updated in the phyloseq object prior to updating its sample data
    new_sam <- sample_data(x) %>% mutate_sample_data(...)
    sample_names(x) <- sample_names(new_sam)
    sample_data(x) <- new_sam
    x
  })

## transmute sample data

#' @rdname mutate-phyloseq
#' @export
setGeneric("transmute_sample_data", 
  function(x, ...) standardGeneric("transmute_sample_data")
)

#' @rdname mutate-phyloseq
setMethod("transmute_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      # Make sure the .sample column is always included; new names can still be
      # set if .sample occurs in the ... expressions
      dplyr::transmute(.sample, ...) %>%
      dplyr::relocate(.sample) %>%
      {suppressMessages(sample_data(.))}
  })

#' @rdname mutate-phyloseq
setMethod("transmute_sample_data", "phyloseq",
  function(x, ...) {
    # Mutating the sample data may create new sample names, which need to get
    # updated in the phyloseq object prior to updating its sample data
    new_sam <- sample_data(x) %>% transmute_sample_data(...)
    sample_names(x) <- sample_names(new_sam)
    sample_data(x) <- new_sam
    x
  })

# rename ----------------------------------------------------------------------

#' Rename columns in the taxonomy table or sample data
#'
#' `rename_tax_table()` and `rename_with_tax_table()` provide the functionality
#' of `dplyr::rename()` and `dplyr::rename_with()` to phyloseq taxonomy tables;
#' `rename_sample_data()` and `rename_with_sample_data()` provide this
#' functionality to phyloseq sample data.
#' 
#' @param x A `phyloseq`, `taxonomyTable`, or `sample_data` object
#' @param ... Renaming expressions passed to `dplyr::rename()` or arguments
#'   passed to `dplyr::rename_with()`
#'
#' @name rename-phyloseq
#'
#' @examples
#' data(GlobalPatterns)
#'
#' GlobalPatterns %>% rank_names
#' 
#' ps1 <- GlobalPatterns %>%
#'   rename_tax_table(Domain = Kingdom) %>%
#'   rename_with_tax_table(stringr::str_to_lower)
#' ps1 %>% rank_names
#' 
#' GlobalPatterns %>% sample_variables
#' ps2 <- GlobalPatterns %>%
#'   rename_with_sample_data(janitor::make_clean_names) %>%
#'   rename_sample_data(sample_id = x_sample_id)
#' ps2 %>% sample_variables
NULL

## rename_tax_table

#' @rdname rename-phyloseq
#' @export
setGeneric("rename_tax_table", 
  function(x, ...) standardGeneric("rename_tax_table")
)

#' @rdname rename-phyloseq
setMethod("rename_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::rename(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname rename-phyloseq
setMethod("rename_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% rename_tax_table(...)
    x
  })

## rename_with_tax_table

#' @rdname rename-phyloseq
#' @export
setGeneric("rename_with_tax_table", 
  function(x, ...) standardGeneric("rename_with_tax_table")
)

#' @rdname rename-phyloseq
setMethod("rename_with_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::rename_with(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname rename-phyloseq
setMethod("rename_with_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% rename_with_tax_table(...)
    x
  })

## rename_sample_data

#' @rdname rename-phyloseq
#' @export
setGeneric("rename_sample_data", 
  function(x, ...) standardGeneric("rename_sample_data")
)

#' @rdname rename-phyloseq
setMethod("rename_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::rename(...) %>%
      {suppressMessages(sample_data(.))}
  })

#' @rdname rename-phyloseq
setMethod("rename_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% rename_sample_data(...)
    x
  })

## rename_with_sample_data

#' @rdname rename-phyloseq
#' @export
setGeneric("rename_with_sample_data", 
  function(x, ...) standardGeneric("rename_with_sample_data")
)

#' @rdname rename-phyloseq
setMethod("rename_with_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::rename_with(...) %>%
      {suppressMessages(sample_data(.))}
  })

#' @rdname rename-phyloseq
setMethod("rename_with_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% rename_with_sample_data(...)
    x
  })
