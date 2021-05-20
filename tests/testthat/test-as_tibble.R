# devtools::test_active_file(here::here("tests", "testthat", "test-as_tibble.R"))

data(GlobalPatterns)
ps <- prune_taxa(
  ape::extract.clade(phy_tree(GlobalPatterns), "0.589.38") %>% .$tip.label,
  GlobalPatterns
)
# add reference sequences
rs <- phy_tree(ps) %>% 
  phangorn::simSeq(10) %>%
  purrr::map(factor, levels = as.character(1:4)) %>%
  purrr::map(forcats::fct_recode, A = "1", C = "2", G = "3", T = "4") %>%
  purrr::map(as.character) %>%
  purrr::map_chr(paste, collapse = "") %>%
  Biostrings::DNAStringSet()
ps <- merge_phyloseq(ps, rs)
rm(rs)


test_that("as_tibble respects name options and handles conflicts", {
  withr::local_options("speedyseq.tibble_sample" = "X.SampleID")
  x <- ps %>% sample_data %>% as_tibble
  expect_identical(names(x)[1:2], c("X.SampleID", "X.SampleID.1"))
  withr::local_options("speedyseq.tibble_otu" = "CL3")
  x <- ps %>% otu_table %>% as_tibble(pivot = FALSE)
  expect_identical(names(x)[1:2], c("CL3", "CL3.1"))
  withr::local_options("speedyseq.tibble_otu" = "X")
  withr::local_options("speedyseq.tibble_abundance" = "X")
  expect_error(x <- ps %>% otu_table %>% as_tibble(pivot = TRUE))
  withr::local_options("speedyseq.tibble_abundance" = "Y")
  x <- ps %>% otu_table %>% as_tibble(pivot = TRUE)
  expect_identical(names(x), c("X", "X.SampleID", "Y"))
  x <- ps %>% as_tibble
  expect_identical(names(x)[1:4], c("X", "X.SampleID", "Y", "X.SampleID.1"))
})
