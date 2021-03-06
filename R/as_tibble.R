#' @importFrom tibble as_tibble
tibble::as_tibble

#' Coerce phyloseq objects to tibble data frames
#'
#' These functions extend the `as_tibble()` function defined by the tibble
#' package to work on phyloseq objects or phyloseq-component objects.
#'
#' Tibbles (tbl_df objects) do not support rownames; the taxa and sample names
#' in the returned tibbles will always be stored in the first or second
#' columns. The names ".otu", ".sample", ".abundance", and ".sequence" are
#' special column names reserved for the otu/taxa names, sample names,
#' abundances, and reference sequences.
#'
#' @param x A phyloseq object or component.
#' @param pivot Whether to pivot the otu table to long format.
#' @param tax Whether to include taxonomy data.
#' @param ref Whether to include reference sequences.
#' @param .name_repair Function to repair names in the case of conflicts.
#'
#' @return A tibble (tbl_df)
#'
#' @name as_tibble-phyloseq
#'
#' @examples
#' library(tibble) # for as_tibble() and glimpse()
#'
#' data(GlobalPatterns)
#'
#' # Subset to 1/100 of the original taxa to speed operations
#' ps <- GlobalPatterns %>%
#'   filter_tax_table(dplyr::row_number() %% 100 == 1)
#'
#' # On phyloseq objects, as_tibble is similar to psmelt()
#' psmelt(ps) %>% glimpse
#' as_tibble(ps) %>% glimpse
#'
#' # By default, the otu_table method provides a tibble in long-format like
#' # psmelt and the phyloseq method
#' otu_table(ps) %>% as_tibble %>% glimpse
#'
#' # Sample data and taxonomy tables produced by as_tibble can be converted
#' # back into their respective phyloseq objects using speedyseq's tbl_df
#' # constructors
#' sample_data(ps) <- sample_data(ps) %>%
#'   as_tibble %>%
#'   dplyr::mutate(sample_sum = sample_sums(ps)) %>%
#'   sample_data
#' sample_data(ps) %>% glimpse
NULL

#' @export
#' @rdname as_tibble-phyloseq
as_tibble.otu_table <- function(x, 
                                pivot = TRUE, 
                                .name_repair = base::make.unique) {
  mat <- x %>% as("matrix")
  if (pivot) {
    if (!taxa_are_rows(x))
      mat <- t(mat)
    tb <- mat %>%
      tibble::as_tibble(rownames = ".otu") %>%
      tidyr::pivot_longer(
        cols = sample_names(x),
        names_to = ".sample",
        values_to = ".abundance"
      )
  } else {
    rn <- ifelse(taxa_are_rows(x), ".otu", ".sample")
    tb <- mat %>%
      tibble::as_tibble(rownames = rn) %>%
      rlang::set_names(., 
        vctrs::vec_as_names(names(.), repair = .name_repair)
      )
  }
  tb
}

#' @export
#' @rdname as_tibble-phyloseq
as_tibble.sample_data <- function(x, .name_repair = base::make.unique) {
  x %>%
    as("data.frame") %>%
    tibble::as_tibble(rownames = ".sample") %>%
    rlang::set_names(., 
      vctrs::vec_as_names(names(.), repair = .name_repair)
    )
}

#' @export
#' @rdname as_tibble-phyloseq
as_tibble.taxonomyTable <- function(x, .name_repair = base::make.unique) {
  x %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = ".otu") %>%
    rlang::set_names(., 
      vctrs::vec_as_names(names(.), repair = .name_repair)
    )
}

#' @export
#' @rdname as_tibble-phyloseq
as_tibble.XStringSet <- function(x) {
  x %>%
    as.character %>%
    tibble::enframe(".otu", ".sequence")
}

#' @export
#' @rdname as_tibble-phyloseq
as_tibble.phyloseq <- function(x, 
                               tax = TRUE, 
                               ref = FALSE, 
                               .name_repair = base::make.unique) {
  tb <- otu_table(x) %>% as_tibble(pivot = TRUE)
  # Add sample data if it exists
  sam <- access(x, "sam_data")
  if (!is.null(sam)) {
    tb <- tb %>% 
      dplyr::left_join(
        sam %>% as_tibble(.name_repair = .name_repair),
        by = ".sample",
        suffix = c("", ".sam")
      )
  }
  # Add tax_table if it exists and tax = TRUE
  tt <- access(x, "tax_table")
  if (tax & !is.null(tt)) {
    tb <- tb %>% 
      dplyr::left_join(
        tt %>% as_tibble(.name_repair = .name_repair),
        by = ".otu", 
        suffix = c("", ".tax")
      )
  }
  # Add refseq if it exists and ref = TRUE
  rs <- access(x, "refseq")
  if (ref & !is.null(rs)) {
    tb <- tb %>% 
      dplyr::left_join(
        rs %>% as_tibble,
        by = ".otu", 
        suffix = c("", ".ref")
      )
  }
  tb
}
