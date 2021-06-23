# devtools::test_active_file(here::here("tests", "testthat", "test-ps_tibble.R"))

data(GlobalPatterns)
ps <- prune_taxa(
  ape::extract.clade(phy_tree(GlobalPatterns), "0.589.38") %>% .$tip.label,
  GlobalPatterns
)
# add reference sequences
rs <- split(
  sample(c("A", "C", "G", "T"), 2 * ntaxa(ps), replace = TRUE),
  rep(seq(ntaxa(ps)), each = 2)
) %>%
  purrr::map_chr(paste, collapse = "") %>%
  rlang::set_names(taxa_names(ps)) %>%
  Biostrings::DNAStringSet()
ps <- merge_phyloseq(ps, rs)
rm(rs)

test_that("ps_tibble returns expected columns", {
  x <- ps %>% sample_data %>% ps_tibble
  expect_identical(names(x), c(".sample", sample_variables(ps)))
  x <- ps %>% otu_table %>% ps_tibble
  expect_identical(names(x), c(".otu", ".sample", ".abundance"))
  x <- ps %>% otu_table %>% ps_tibble(pivot = FALSE)
  expect_identical(names(x), c(".otu", sample_names(ps)))
  x <- ps %>% ps_tibble(ref = TRUE)
  expect_identical(names(x), 
    c(".otu", ".sample", ".abundance", sample_variables(ps), rank_names(ps),
      ".sequence")
  )
})

test_that("ps_tibble and psmelt agree up to row order", {
  x <- ps %>% ps_tibble %>% 
    dplyr::rename(OTU = .otu, Sample = .sample, Abundance = .abundance)
  y <- ps %>% psmelt
  expect_true(dplyr::all_equal(x, y, ignore_row_order = TRUE))
})
