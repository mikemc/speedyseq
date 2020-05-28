# run just this file:
# devtools::test_file(here::here("tests", "testthat", "test-select.R"))

library(magrittr)

# Get a phyloseq object with all slots filled for testing
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
# Taxa to select, in order different from original
taxa <- taxa_names(ps)[c(5,1,9,2,6,3)]
current_order <- taxa_names(ps) %>% intersect(taxa)

test_that("select_taxa() gives taxa in requested order unless `x` is or contains a tree", {
  expect_equal(current_order, ps %>% select_taxa(taxa) %>% taxa_names)
  expect_equal(current_order, ps %>% phy_tree %>% select_taxa(taxa) %>% taxa_names)
  expect_equal(taxa, ps %>% refseq %>% select_taxa(taxa) %>% taxa_names)
  expect_equal(taxa, ps %>% otu_table %>% select_taxa(taxa) %>% taxa_names)
  expect_equal(taxa, ps %>% tax_table %>% select_taxa(taxa) %>% taxa_names)
  expect_equal(current_order, 
    ps %>% refseq %>% select_taxa(taxa, reorder = FALSE) %>% taxa_names)
  expect_equal(current_order, 
    ps %>% otu_table %>% select_taxa(taxa, reorder = FALSE) %>% taxa_names)
  expect_equal(current_order, 
    ps %>% tax_table %>% select_taxa(taxa, reorder = FALSE) %>% taxa_names)

  ps1 <- ps
  ps1@phy_tree <- NULL
  expect_equal(taxa, ps1 %>% select_taxa(taxa) %>% taxa_names)
  expect_equal(current_order, 
    ps1 %>% select_taxa(taxa, reorder = FALSE) %>% taxa_names)
  ps1@refseq <- NULL
  expect_equal(taxa, ps1 %>% select_taxa(taxa) %>% taxa_names)
  ps1@tax_table <- NULL
  expect_equal(taxa, ps1 %>% select_taxa(taxa) %>% taxa_names)
})

test_that("select_taxa() agrees with prune_taxa() when accounting for order", {
  expect_equal(
    ps %>% select_taxa(taxa) %>% otu_table,
    ps %>% prune_taxa(taxa, .) %>% otu_table
  )
  expect_equal(
    ps %>% sample_data %>% select_taxa(taxa),
    ps %>% sample_data %>% prune_taxa(taxa, .)
  )
  expect_equal(
    ps %>% phy_tree %>% select_taxa(taxa),
    ps %>% phy_tree %>% prune_taxa(taxa, .)
  )
  expect_equal(
    ps %>% refseq %>% select_taxa(taxa, reorder = FALSE),
    ps %>% refseq %>% prune_taxa(taxa, .)
  )
  expect_equal(
    ps %>% otu_table %>% select_taxa(taxa, reorder = TRUE),
    ps %>% otu_table %>% prune_taxa(taxa, .)
  )
  expect_equal(
    ps %>% tax_table %>% select_taxa(taxa, reorder = TRUE),
    ps %>% tax_table %>% prune_taxa(taxa, .)
  )
})

test_that("select_taxa() returns sample data unchanged", {
  expect_equal(ps %>% sample_data %>% select_taxa(taxa), sample_data(ps))
})

test_that("select_taxa() throws an error on duplicates or taxa not in x", {
  expect_error(ps %>% select_taxa(c(taxa, taxa[1])),
    "is not TRUE"
  )
  expect_error(ps %>% sample_data %>% select_taxa(c(taxa, taxa[1])),
    "is not TRUE"
  )
  expect_error(ps %>% otu_table %>% select_taxa(c(taxa, taxa[1])),
    "is not TRUE"
  )
  expect_error(ps %>% phy_tree %>% select_taxa(c(taxa, taxa[1])),
    "is not TRUE"
  )
  expect_error(ps %>% tax_table %>% select_taxa(c(taxa, taxa[1])),
    "is not TRUE"
  )
  expect_error(ps %>% refseq %>% select_taxa(c(taxa, taxa[1])),
    "is not TRUE"
  )
  expect_error(ps %>% select_taxa(c(taxa, "asdf")),
    "is not TRUE"
  )
  expect_error(ps %>% otu_table %>% select_taxa(c(taxa, "asdf")),
    "is not TRUE"
  )
  expect_error(ps %>% phy_tree %>% select_taxa(c(taxa, "asdf")),
    "is not TRUE"
  )
  expect_error(ps %>% tax_table %>% select_taxa(c(taxa, "asdf")),
    "is not TRUE"
  )
  expect_error(ps %>% refseq %>% select_taxa(c(taxa, "asdf")),
    "is not TRUE"
  )
})

