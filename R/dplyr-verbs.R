# vim: foldmethod=marker

# Some of the documentation for these functions is modified from the
# corresponding dplyr functions (https://github.com/tidyverse/dplyr), MIT
# license RStudio and others.

# filter -------------------------------------------------------------------{{{

#' Subset rows in the taxonomy table or sample data using column values
#'
#' These functions are wrappers around `dplyr::filter()` that make it possible
#' to subset rows (corresponding to taxa or samples) using column values in the
#' taxonomy table or sample data.
#' 
#' @param x A `phyloseq`, `taxonomyTable`, or `sample_data` object
#' @param ... Expressions passed to `dplyr::filter()`
#'
#' @name filter-phyloseq
#'
#' @examples
#' data(GlobalPatterns)
#'
#' ps <- GlobalPatterns %>%
#'   filter_tax_table(Kingdom == "Bacteria") %>%
#'   filter_sample_data(SampleType %in% c("Feces", "Soil"))
NULL

#' @rdname filter-phyloseq
#' @export
setGeneric("filter_tax_table", 
  function(x, ...) standardGeneric("filter_tax_table")
)

#' @rdname filter-phyloseq
setMethod("filter_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% filter_tax_table(...)
    x
  })

#' @rdname filter-phyloseq
setMethod("filter_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::filter(...) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname filter-phyloseq
#' @export
setGeneric("filter_sample_data", 
  function(x, ...) standardGeneric("filter_sample_data")
)

#' @rdname filter-phyloseq
setMethod("filter_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% filter_sample_data(...)
    x
  })

#' @rdname filter-phyloseq
setMethod("filter_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::filter(...) %>%
      {suppressMessages(sample_data(.))}
  })

# }}}

# mutate and transmute -----------------------------------------------------{{{

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

# }}}

# select -------------------------------------------------------------------{{{

#' Subset columns in the taxonomy table or sample data using their names and types
#'
#' These functions are wrappers around `dplyr::select()` that provide
#' convenient ways to modify `tax_table(x)` and `sample_data(x)`.
#' See `dplyr::select()` for supported syntax and helpers.
#'
#' The special column names '.otu' and '.sample' should not be used; see
#' `mutate-phyloseq` for the ability to change taxa and sample names using
#' these names.
#'
#' @param x A `phyloseq`, `taxonomyTable`, or `sample_data` object
#' @param ... Expressions passed to `dplyr::select()`
#'
#' @name select-phyloseq
#'
#' @examples
#' data(GlobalPatterns)
#'
#' GlobalPatterns %>% rank_names
#' GlobalPatterns %>% sample_variables
#' ps <- GlobalPatterns %>%
#'   select_tax_table(Phylum, Genus:Species) %>%
#'   select_sample_data(!dplyr::contains("Barcode"))
#' ps %>% rank_names
#' ps %>% sample_variables
NULL

#' @rdname select-phyloseq
#' @export
setGeneric("select_tax_table", 
  function(x, ...) standardGeneric("select_tax_table")
)

#' @rdname select-phyloseq
setMethod("select_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      # Ensure .otu always kept; Putafter the ... to handle use of '-' for
      # column removal
      dplyr::select(..., .otu) %>%
      dplyr::relocate(.otu) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname select-phyloseq
setMethod("select_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% select_tax_table(...)
    x
  })

#' @rdname select-phyloseq
#' @export
setGeneric("select_sample_data", 
  function(x, ...) standardGeneric("select_sample_data")
)

#' @rdname select-phyloseq
setMethod("select_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      # Ensure .sample always kept; put after the ... to handle use of '-' for
      # column removal
      dplyr::select(..., .sample) %>%
      dplyr::relocate(.sample) %>%
      {suppressMessages(sample_data(.))}
  })

#' @rdname select-phyloseq
setMethod("select_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% select_sample_data(...)
    x
  })

# }}}

# relocate -----------------------------------------------------------------{{{

#' Change column order in the taxonomy table or sample data
#'
#' These functions are wrappers around `dplyr::relocate()` that provide
#' convenient ways to modify `tax_table(x)` and `sample_data(x)`.
#'
#' The special column names '.otu' and '.sample' should not be used; see
#' `mutate-phyloseq` for the ability to change taxa and sample names using
#' these names.
#'
#' @param x A `phyloseq`, `taxonomyTable`, or `sample_data` object
#' @param ... Expressions and arguments passed to `dplyr::relocate()`
#'
#' @name relocate-phyloseq
#'
#' @examples
#' data(GlobalPatterns)
#'
#' GlobalPatterns %>% sample_variables
#' ps <- GlobalPatterns %>%
#'   relocate_sample_data(SampleType, .after = X.SampleID)
#' ps %>% sample_variables
NULL

#' @rdname relocate-phyloseq
#' @export
setGeneric("relocate_tax_table", 
  function(x, ...) standardGeneric("relocate_tax_table")
)

#' @rdname relocate-phyloseq
setMethod("relocate_tax_table", "taxonomyTable",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::relocate(...) %>%
      # Ensure .otu always kept first
      dplyr::relocate(.otu) %>%
      {suppressMessages(tax_table(.))}
  })

#' @rdname relocate-phyloseq
setMethod("relocate_tax_table", "phyloseq",
  function(x, ...) {
    tax_table(x) <- tax_table(x) %>% relocate_tax_table(...)
    x
  })

#' @rdname relocate-phyloseq
#' @export
setGeneric("relocate_sample_data", 
  function(x, ...) standardGeneric("relocate_sample_data")
)

#' @rdname relocate-phyloseq
setMethod("relocate_sample_data", "sample_data",
  function(x, ...) {
    x %>%
      ps_tibble %>%
      dplyr::relocate(...) %>%
      # Ensure .sample always kept first
      dplyr::relocate(.sample) %>%
      {suppressMessages(sample_data(.))}
  })

#' @rdname relocate-phyloseq
setMethod("relocate_sample_data", "phyloseq",
  function(x, ...) {
    sample_data(x) <- sample_data(x) %>% relocate_sample_data(...)
    x
  })

# }}}

# rename -------------------------------------------------------------------{{{

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
