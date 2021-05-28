#' Coerce phyloseq objects to tibbles (tbl_df objects)
#'
#' An S4 generic function that serves the role of `tibble::as_tibble()` on
#' phyloseq and phyloseq-component objects.
#'
#' @param x A phyloseq object or component.
#' @param pivot Whether to pivot the data frame to long format.
#' @param tax Whether to include taxonomy data.
#' @param ref Whether to include reference sequences.
#' @param .name_repair Function to repair names in the case of conflicts.
#' @param ... Unused; for methods.
#'
#' @name ps_tibble
setGeneric("ps_tibble", 
  function(x, ...) standardGeneric("ps_tibble")
)

#' @rdname ps_tibble
setMethod("ps_tibble", "otu_table",
  function(x, pivot = TRUE, .name_repair = base::make.unique) {
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
  })

#' @rdname ps_tibble
setMethod("ps_tibble", "sample_data",
  function(x, .name_repair = base::make.unique) {
    x %>%
      as("data.frame") %>%
      tibble::as_tibble(rownames = getOption("speedyseq.tibble_sample")) %>%
      rlang::set_names(., 
        vctrs::vec_as_names(names(.), repair = .name_repair)
      )
  })

#' @rdname ps_tibble
setMethod("ps_tibble", "taxonomyTable",
  function(x, .name_repair = base::make.unique) {
    x %>%
      as("matrix") %>%
      tibble::as_tibble(rownames = getOption("speedyseq.tibble_otu")) %>%
      rlang::set_names(., 
        vctrs::vec_as_names(names(.), repair = .name_repair)
      )
  })

#' @rdname ps_tibble
setMethod("ps_tibble", "XStringSet",
  function(x) {
    nms <- c(
      getOption("speedyseq.tibble_otu"), 
      getOption("speedyseq.tibble_sequence")
    )
    if (anyDuplicated(nms)) 
      stop("OTU and sequence column names must be unique")
    x %>%
      as.character %>%
      tibble::enframe(nms[1], nms[2])
  })

#' @rdname ps_tibble
setMethod("ps_tibble", "phyloseq",
  function(x, tax = TRUE, ref = FALSE, .name_repair = base::make.unique) {
    nms <- c(
      getOption("speedyseq.tibble_otu"), 
      getOption("speedyseq.tibble_sample"), 
      getOption("speedyseq.tibble_abundance")
    )
    if (anyDuplicated(nms)) 
      stop("OTU, sample, and abundance column names must be unique")
    tb <- otu_table(x) %>% ps_tibble(pivot = TRUE)
    # Add sample data if it exists
    sam <- access(x, "sam_data")
    if (!is.null(sam)) {
      tb <- tb %>% 
        dplyr::left_join(
          sam %>% ps_tibble(.name_repair = .name_repair),
          by = nms[2],
          suffix = c("", ".sam")
        )
    }
    # Add tax_table if it exists and tax = TRUE
    tt <- access(x, "tax_table")
    if (tax & !is.null(tt)) {
      tb <- tb %>% 
        dplyr::left_join(
          tt %>% ps_tibble(.name_repair = .name_repair),
          by = nms[1], 
          suffix = c("", ".tax")
        )
    }
    # Add refseq if it exists and ref = TRUE
    rs <- access(x, "refseq")
    if (ref & !is.null(rs)) {
      tb <- tb %>% 
        dplyr::left_join(
          rs %>% ps_tibble,
          by = nms[1], 
          suffix = c("", ".ref")
        )
    }
    tb
  })
