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

test_that("ps_tibble respects name options and handles conflicts", {
  withr::local_options("speedyseq.tibble_sample" = "X.SampleID")
  x <- ps %>% sample_data %>% ps_tibble
  expect_identical(names(x)[1:2], c("X.SampleID", "X.SampleID.1"))
  withr::local_options("speedyseq.tibble_otu" = "CL3")
  x <- ps %>% otu_table %>% ps_tibble(pivot = FALSE)
  expect_identical(names(x)[1:2], c("CL3", "CL3.1"))
  withr::local_options("speedyseq.tibble_otu" = "X")
  withr::local_options("speedyseq.tibble_abundance" = "X")
  expect_error(x <- ps %>% otu_table %>% ps_tibble(pivot = TRUE))
  withr::local_options("speedyseq.tibble_abundance" = "Y")
  x <- ps %>% otu_table %>% ps_tibble(pivot = TRUE)
  expect_identical(names(x), c("X", "X.SampleID", "Y"))
  x <- ps %>% ps_tibble
  expect_identical(names(x)[1:4], c("X", "X.SampleID", "Y", "X.SampleID.1"))
})
