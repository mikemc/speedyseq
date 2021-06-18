#' Merge samples by a sample variable or factor
#'
#' This function provides an alternative to `phyloseq::merge_samples()` that
#' better handles sample variables of different types, especially categorical
#' sample variables. It combines the samples in `x` defined by the sample
#' variable or factor `group` by summing the abundances in `otu_table(x)` and
#' combines sample variables by the summary functions in `funs`. The default
#' summary function, `unique_or_na()`, collapses the values within a group to a
#' single unique value if it exists and otherwise returns NA. The new (merged)
#' samples are named by the values in `group`.
#' 
#' @param x A `phyloseq`, `otu_table`, or `sample_data` object
#' @param group A sample variable or a vector of length `nsamples(x)` defining
#'   the sample grouping. A vector must be supplied if x is an otu_table
#' @param fun_otu Function for combining abundances in the otu table; default
#'   is `sum`. Can be a formula to be converted to a function by
#'   [purrr::as_mapper()]
#' @param funs Named list of merge functions for sample variables; default is
#'   `unique_or_na`
#' @param reorder Logical specifying whether to reorder the new (merged)
#'   samples by name
#' 
#' @export
#'
#' @examples
#' data(enterotype)
#' 
#' # Merge samples with the same project and clinical status
#' ps <- enterotype 
#' sample_data(ps) <- sample_data(ps) %>%
#'   transform(Project.ClinicalStatus = Project:ClinicalStatus)
#' sample_data(ps) %>% head
#' ps0 <- merge_samples2(ps, "Project.ClinicalStatus",
#'   fun_otu = mean,
#'   funs = list(Age = mean)
#' )
#' sample_data(ps0) %>% head
setGeneric("merge_samples2", 
  function(x,
           group,
           fun_otu = sum,
           funs = list(),
           reorder = FALSE)
    standardGeneric("merge_samples2")
)

#' @rdname merge_samples2
setMethod(
  "merge_samples2", 
  signature("phyloseq"), 
  function(x, group, fun_otu = sum, funs = list(), reorder = FALSE) {
    if (length(group) == 1) {
      stopifnot(group %in% sample_variables(x))
      group <- sample_data(x)[[group]]
    } else {
      stopifnot(identical(length(group), nsamples(x)))
    }
    # Drop samples with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding samples will be dropped")
      x <- prune_samples(!is.na(group), x)
      group <- group[!is.na(group)]
    }
    # Merge
    otu.merged <- merge_samples2(otu_table(x), group, 
      fun_otu = fun_otu, 
      reorder = reorder
    )
    if (!is.null(access(x, "sam_data")))
      sam.merged <- merge_samples2(sample_data(x), group, funs = funs)
    else 
      sam.merged <- NULL
    phyloseq(
      otu.merged, 
      sam.merged,
      access(x, "tax_table"),
      access(x, "phy_tree"),
      access(x, "refseq")
    )
})

#' @rdname merge_samples2
setMethod(
  "merge_samples2", 
  signature("otu_table"), 
  function(x, group, fun_otu = sum, reorder = FALSE) {
    stopifnot(identical(length(group), nsamples(x)))
    # Work with samples as rows, and remember to flip back at end if needed
    needs_flip <- taxa_are_rows(x)
    if (needs_flip)
      x <- t(x)
    # Drop samples with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding samples will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    # Merging; result is a matrix with taxa as columns and rownames
    # corresponding to `group`
    if (identical(fun_otu, sum)) {
      x.merged <- rowsum(x, group, reorder = reorder)
    } else {
      stopifnot(!".group" %in% colnames(x))
      f <- purrr::as_mapper(fun_otu)
      x <- x %>% 
        as("matrix") %>%
        data.table::as.data.table() %>%
        cbind(.group = group)
      if (reorder) 
        x.merged <- x[, lapply(.SD, f), keyby = .(.group)]
      else
        x.merged <- x[, lapply(.SD, f), by = .(.group)]
      rns <- x.merged$.group
      x.merged[, .group := NULL]
      x.merged <- x.merged %>% as("matrix")
      rownames(x.merged) <- rns
    }
    # Return an otu table in the proper orientation
    x.merged <- x.merged %>% otu_table(taxa_are_rows = FALSE)
    if (needs_flip)
      x.merged <- t(x.merged)
    x.merged
})

#' @rdname merge_samples2
setMethod(
  "merge_samples2", 
  signature("sample_data"), 
  function(x, group, funs = list(), reorder = FALSE) {
    if (length(group) == 1) {
      stopifnot(group %in% sample_variables(x))
      group <- x[[group]]
    } else {
      stopifnot(identical(length(group), nsamples(x)))
    }
    # Drop samples with `is.na(group)`
    if (anyNA(group)) {
      warning("`group` has missing values; corresponding samples will be dropped")
      x <- x[!is.na(group), ]
      group <- group[!is.na(group)]
    }
    ## Set the functions f used to merge each sample variable.
    # Named logical vector indicating whether each variable is in the funs
    var_in_funs <- names(x) %>% 
      rlang::set_names(. %in% names(funs), .)
    # For vars in the funs, run f through as_mapper; else, use the default f
    funs <- purrr::map2(var_in_funs, names(var_in_funs),
      ~if (.x) purrr::as_mapper(funs[[.y]]) else unique_or_na
    )
    ## Merge variable values, creating a new sample_data object with one row
    ## per group.
    # A "sample_data" object is a list of data variables (columns); strategy is
    # to reduce each variable with `merge_groups()`, and then recombine into a
    # data.frame. The call to `merge_groups()` will sort by `group` values,
    # which we need to account for when setting the new sample names.
    new_sample_names <- group %>% unique %>% sort %>% as.character
    x.merged <- purrr::map2(x, funs, 
      ~merge_groups(.x, group = group, f = .y)
    ) %>%
      data.frame %>%
      vctrs::vec_set_names(new_sample_names)
    ## Put back in initial order
    if (!reorder) {
      initial_order <- group %>% unique %>% as.character
      x.merged <- x.merged[initial_order, , drop = FALSE]
    }
    ## Return as sample data with group names preserved
    x.merged %>% sample_data_stable
})

# Helpers ---------------------------------------------------------------------

#' Get the unique value in x or NA if none
#'
#' If `unique(x)` is a single value, return it; otherwise, return an NA of the
#' same type as `x`. If `x` is a factor, then the levels and ordered status
#' will be kept in either case. If `x` is a non-atomic vector (i.e. a list),
#' then the logical `NA` will be used.
#'
#' @param x A vector
#' @export
#' 
#' @examples
#' f <- factor(c("a", "a", "b", "c"), ordered = TRUE)
#' unique_or_na(f)
#' unique_or_na(f[1:2])
#' 
#' x <- c("a", "b", "a")
#' unique_or_na(x[c(1, 3)])
#' unique_or_na(x)
#' unique_or_na(x) %>% typeof
unique_or_na <- function(x) {
  UseMethod("unique_or_na")
}

#' @export
unique_or_na.default <- function(x) {
  if (length(unique(x)) == 1)
    x[[1]]
  else if (is.atomic(x))
    as(NA, typeof(x))
  else
    NA
}

#' @export
unique_or_na.factor <- function(x) {
  if (length(unique(x)) == 1)
    x[[1]]
  else
    factor(NA, levels = levels(x), ordered = is.ordered(x))
}

#' Merge groups of elements within a vector by a function
#'
#' Internal function used in `merge_samples2()` to merge variables. Note, owing
#' to the use of `split()`, the merged elements in the new vector will be
#' reordered according to `group`.
#'
#' @param x A vector whose elements will be merged.
#' @param group A vector such that `as.factor(group)` defines the grouping.
#' @param f A function that, when applied to a subvector of x, returns a single
#'   value. Can also be a formula as interpretted by `purrr::as_mapper()`.
merge_groups <- function(x, group, f = unique_or_na) {
  f <- purrr::as_mapper(f)
  split(x, group) %>% 
    purrr::map(f) %>%
    {vctrs::vec_c(!!!., .name_spec = rlang::zap())}
}

#' Create sample data without adjusting row/sample names
#'
#' `phyloseq::sample_data()` will change the sample names from the row names if
#' they are `as.character(1:nrow(object))`. This function instead keeps the
#' names as is.
#' 
#' @param object A "data.frame"-class object
#'
#' @keywords internal
#'
#' @examples 
#' x <- data.frame(var1 = letters[1:3], var2 = 7:9)
#' rownames(x)
#' sample_data(x)
#' speedyseq:::sample_data_stable(x)
sample_data_stable <- function(object) {
  # Modified from phyloseq's sample_data data.frame method; see
  # https://github.com/joey711/phyloseq/blob/master/R/sampleData-class.R
  stopifnot(identical(class(object), "data.frame"))
	# Make sure there are no phantom levels in categorical variables
	object <- phyloseq:::reconcile_categories(object)
	# instantiate first to check validity
	SM <- new("sample_data", object)
	SM
}
