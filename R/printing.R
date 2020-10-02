# show/print ------------------------------------------------------------------

# Phyloseq objects

#' Methods for `show()` for phyloseq objects.
#'
#' See \code{\link[methods]{show}}. Speedyseq provides `show()` methods that
#' overwrite phyloseq's methods for the classes defined in phyloseq.
#'
#' @seealso \code{\link[methods]{show}}
#' 
#' @inheritParams methods::show
#' @examples
#' data(GlobalPatterns)
#' show(GlobalPatterns)
#' GlobalPatterns
#' @rdname show-methods
#' @export
setMethod("show", "phyloseq", 
  function(object){
    cat("phyloseq-class experiment-level object", fill = TRUE)
    # print otu_table (always there).
    paste0("otu_table()   ", otu_table(object) %>% oneline) %>%
      cat(fill = TRUE)
    # print Sample Data if there
    if(!is.null(sample_data(object, FALSE))){
      paste0("sample_data() ", sample_data(object) %>% oneline) %>%
        cat(fill = TRUE)
    }
    # print tax table if there	
    if(!is.null(tax_table(object, FALSE))){
      paste0("tax_table()   ", tax_table(object) %>% oneline) %>%
        cat(fill = TRUE)
    }
    # print tree if there
    if(!is.null(phy_tree(object, FALSE))){
      paste0("phy_tree()    ", phy_tree(object) %>% oneline) %>%
        cat(fill = TRUE)
    }
    # print refseq summary if there
    if(!is.null(refseq(object, FALSE))){
      paste0("refseq()      ", refseq(object) %>% oneline) %>%
        cat(fill = TRUE)
    }
    if (taxa_are_rows(object))
      cat("taxa are rows", fill = TRUE)
    else
      cat("taxa are columns", fill = TRUE)
  }
)

# Sample data and tax table, using tibble_print()

#' @rdname show-methods
setMethod("show", "taxonomyTable",
  function(object) {
    cat(oneline(object), fill = TRUE)
    object %>%
      as("data.frame") %>%
      tibble_print(rows = "taxa", cols = "taxonomic ranks")
  }
)

#' @export
print.taxonomyTable <- function(x, 
                                ..., 
                                n = NULL, 
                                width = NULL, 
                                n_extra = NULL) {
  cat(oneline(x), fill = TRUE)
  x %>%
    as("matrix") %>%
    tibble_print(n = n, width = width, n_extra = n_extra, types = FALSE,
      rows = "taxa", cols = "taxonomic ranks"
    )
  invisible(x)
}

#' @rdname show-methods
setMethod("show", "sample_data", 
  function(object) {
    cat(oneline(object), fill = TRUE)
    object %>%
      as("data.frame") %>%
      tibble_print(types = FALSE, rows = "samples", cols = "variables")
  }
)

#' @export
print.sample_data <- function(x, 
                              ..., 
                              n = NULL, 
                              width = NULL, 
                              n_extra = NULL) {
  cat(oneline(x), fill = TRUE)	
  x %>%
    as("data.frame") %>%
    tibble_print(n = n, width = width, n_extra = n_extra,
      rows = "samples", cols = "variables"
    )
  invisible(x)
}

#' @rdname show-methods
setMethod("show", "otu_table", 
  function(object) {
    cat(oneline(object), fill = TRUE)	
    if (taxa_are_rows(object)) {
      cat("Taxa are rows", fill = TRUE)
      rows <- "taxa (rows)"
      cols <- "samples (columns)"
    }
    else {
      cat("Taxa are columns", fill = TRUE)
      cols <- "taxa (columns)"
      rows <- "samples (rows)"
    }
    object %>%
      as("matrix") %>%
      tibble_print(n_extra = 0, types = FALSE, rows = rows, cols = cols)
  }
)

#' @export
print.otu_table <- function(x, 
                            ..., 
                            n = NULL, 
                            width = NULL, 
                            n_extra = 0) {
  cat(oneline(x), fill = TRUE)
  if (taxa_are_rows(x)) {
    cat("Taxa are rows", fill = TRUE)
    rows <- "taxa (rows)"
    cols <- "samples (columns)"
  }
  else {
    cat("Taxa are columns", fill = TRUE)
    cols <- "taxa (columns)"
    rows <- "samples (rows)"
  }
  x %>%
    as("matrix") %>%
    tibble_print(n = n, width = width, n_extra = n_extra, types = FALSE,
      rows = rows, cols = cols)
  invisible(x)
}

# Helpers ---------------------------------------------------------------------

setAs("taxonomyTable", "data.frame",
  function(from) {
    from %>% as("matrix") %>% data.frame
  }
)

#' Print a 2-d data.frame or matrix like a tibble
#'
#' @param x A data frame or matrix
#' @param n Number of rows to print
#' @param width Display width
#' @param n_extra Number of extra columns after width to summarize after main
#'   display
#' @param types Whether to print the second row with the variable types
#' @param rows What rows are called in the end summary
#' @param cols What columns are called in the end summary
#'
#' @keywords internal
tibble_print <- function(x, 
                         n = NULL, 
                         width = NULL, 
                         n_extra = NULL, 
                         types = TRUE,
                         rows = "rows", 
                         cols = "variables") {
  stopifnot(is.data.frame(x) || is.matrix(x))
  # Set n; Modified from tibble::trunc_mat
  nrows <- nrow(x)
  if (is.null(n) || n < 0) {
    if (is.na(nrows) || nrows > 20) {
      n <- 10
    } else {
      n <- nrows
    }
  }
  # How many characters will be taken up by the row numbers?
  nchar_num <- nchar(as.character(n)) + 1
  if (is.null(width))
    width <- getOption("width")
  # Add a few chars to account for the row numbers that will be cut
  width <- width + nchar_num
  # Also increase the `width` setting to avoid wrapping
  rlang::scoped_options(width = getOption("width") + nchar_num)
  # the first n rownames will be printed in place of the row numbers
  rns <- rownames(x) %>% utils::head(n)
  # how many characters to show of the rownames column. 
  nchar_rn <- min(max(nchar(rns), nchar("<chr>")), 10)
  # column name for the rownames
  rn <- stringr::str_pad("rn", nchar_rn, "left", ".")
  # Print x as a tibble, and save the output in s. Limit number of columns of x
  # to speed up call to trunc_mat when ncol(x) is very large
  ncol_x <- ncol(x)
  x <- x[, seq(min(ncol_x, as.integer(width/5 + n_extra + 1))), drop = FALSE]
  s <- x %>%
    tibble::as_tibble(rownames = rn) %>%
    {utils::capture.output(
      tibble::trunc_mat(., n = n, width = width, n_extra = n_extra)
    )} %>%
    utils::tail(-1)
  # Split off last section with the additional rows and columns
  last_section <- s[stringr::str_detect(s, "^#")]
  s <- s[!stringr::str_detect(s, "^#")]
  # Adjust the last section to account for any cut columns and to use the
  # correct `rows` and `cols` names
  ncol_more <- last_section %>% 
    stringr::str_extract("(?<=, and )[0-9]+(?= more variables)") %>%
    as.integer
  ncol_shown <- ncol(x) - ncol_more
  ncol_more_actual <- as.character(ncol_x - ncol_shown)
  last_section <- last_section %>%
    stringr::str_replace("(?<=, and )[0-9]+(?= more variables)", 
      ncol_more_actual) %>%
    stringr::str_replace("rows", rows) %>%
    stringr::str_replace("variables", cols)
  # Remove the row numbers. To keep, comment out and add nchar_num nchar_rn + 1
  # in the next two str_sub commands
  s <- purrr::map_chr(s, stringr::str_sub, nchar_num + 1)
  # Remove the col name and type in the first col
  w <- nchar(s[1])
  s[1] <- stringr::str_sub(s[1], nchar_rn + 1) %>% stringr::str_pad(w, "left")
  s[2] <- stringr::str_sub(s[2], nchar_rn + 1) %>% stringr::str_pad(w, "left")
  if (!types)
    s <- s[-2]
  cat(s, sep = "\n")
  cat(last_section, sep = "\n")
}


#' One-line summaries of phyloseq components
#'
#' @param x A phyloseq component object
#'
#' @keywords internal
oneline <- function(x) {
  UseMethod("oneline")
}

oneline.otu_table <- function(x) {
  stringr::str_glue(
    "OTU Table:          [ {ntaxa(x)} taxa and {nsamples(x)} samples ]:", 
  )
}

oneline.sample_data <- function(x) {
  stringr::str_glue(
    "Sample Data:        [ {dim(x)[1]} samples by {dim(x)[2]} sample variables ]:", 
  )
}

oneline.taxonomyTable <- function(x) {
  stringr::str_glue(
    "Taxonomy Table:     [ {dim(x)[1]} taxa by {dim(x)[2]} taxonomic ranks ]:", 
  )
}

oneline.phylo <- function(x) {
  stringr::str_glue(
    "Phylogenetic Tree:  [ {ntaxa(x)} tips and {phy_tree(x)$Nnode} internal nodes ]:", 
  )
}

oneline.XStringSet <- function(x) {
  stringr::str_glue(
    "{class(x)}:       [ {ntaxa(x)} reference sequences ]"
  )
}

setMethod("oneline", "otu_table", oneline.otu_table)
setMethod("oneline", "sample_data", oneline.sample_data)
setMethod("oneline", "taxonomyTable", oneline.taxonomyTable)
setMethod("oneline", "phylo", oneline.phylo)
setMethod("oneline", "XStringSet", oneline.XStringSet)

# Methods for dplyr::glimpse --------------------------------------------------

# Fix output of tibble::glimpse() on sample_data and taxonomyTable objects by
# first converting to a data.frame. Note, currently the sample and taxa names
# will not be shown.

#' @importFrom tibble glimpse
#' @export
glimpse.sample_data <- function(x, width = NULL, ...) {
  x %>% 
    as("data.frame") %>% 
    tibble::glimpse(width = width, ...)
}

#' @importFrom tibble glimpse
#' @export
glimpse.taxonomyTable <- function(x, width = NULL, ...) {
  x %>% 
    as.data.frame %>% 
    tibble::glimpse(width = width, ...)
}


