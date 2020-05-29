# run just this file:
# devtools::test_file(here::here("tests", "testthat", "test-glom-functions.R"))
context("Equivalence of tax_glom() with phyloseq::tax_glom()")

library(magrittr) # Imported for the pipe

data(GlobalPatterns)

# tax_glom--------------------------------------------------------------------- 

# Pick a subset to speed up the tests
set.seed(20190421)
ps <- GlobalPatterns %>%
    {prune_taxa(sample(taxa_names(.), 2e2), .)}
tax <- tax_table(ps)
# NOTE: taxon order is fixed by the tree in this dataset.

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for an intermediate rank with default arguments", {
  ps1 <- phyloseq::tax_glom(ps, "Phylum")
  ps2 <- tax_glom(ps, "Phylum")
  expect_identical(ps1, ps2)
  # Also test on a tax_table only
  tax1 <- phyloseq::tax_glom(tax, "Phylum")
  tax2 <- tax_glom(tax, "Phylum")
  expect_identical(tax1, tax2)
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for an intermediate rank with NArm = FALSE", {
  ps1 <- phyloseq::tax_glom(ps, "Family", NArm = FALSE)
  ps2 <- tax_glom(ps, "Family", NArm = FALSE)
  expect_identical(ps1, ps2)
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for the highest rank", {
  ps1 <- phyloseq::tax_glom(ps, "Kingdom")
  ps2 <- tax_glom(ps, "Kingdom")
  expect_identical(ps1, ps2)
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for the lowest rank", {
  ps1 <- phyloseq::tax_glom(ps, "Species")
  ps2 <- tax_glom(ps, "Species")
  expect_identical(ps1, ps2)
})

test_that("tax_glom() gives taxa in expected order", {
  idx <- c(5000, 18000, 1, 10000, 2, 5001)
  taxa <- taxa_names(GlobalPatterns)[idx]
  ps0 <- select_taxa(GlobalPatterns, taxa)
  # Reorder taxa after dropping tree
  ps0@phy_tree <- NULL
  ps0 <- select_taxa(ps0, taxa)
  # tibble for getting the expected results
  tb <- tax_table(ps0) %>%
    as("matrix") %>%
    tibble::as_tibble(rownames = "otu") %>%
    tibble::add_column(sum = taxa_sums(ps0)) %>%
    dplyr::mutate(row = dplyr::row_number()) %>%
    dplyr::group_by(Kingdom, Phylum) %>%
    dplyr::summarize(
      max_otu = otu[which.max(sum)],
      first_otu = otu[which.min(row)],
      first_row = min(row)
    ) %>%
    dplyr::ungroup() %>%
    tibble::add_column(.before = 3,
      Class = NA_character_, 
      Order = NA_character_, 
      Family = NA_character_, 
      Genus = NA_character_, 
      Species = NA_character_
    )
  ## tax_glom on phyloseq object, reorder = FALSE
  ps1 <- tax_glom(ps0, "Phylum", reorder = FALSE)
  # expected tax table
  tt1 <- tb %>%
    dplyr::arrange(first_row) %>%
    dplyr::select(max_otu, rank_names(ps0)) %>%
    tibble::column_to_rownames("max_otu") %>%
    as("matrix") %>%
    tax_table
  expect_identical(tax_table(ps1), tt1)
  ## tax_glom on phyloseq object, reorder = TRUE
  ps2 <- tax_glom(ps0, "Phylum", reorder = TRUE)
  # expected tax table
  tt2 <- tb %>%
    dplyr::arrange(Kingdom, Phylum) %>%
    dplyr::select(max_otu, rank_names(ps0)) %>%
    tibble::column_to_rownames("max_otu") %>%
    as("matrix") %>%
    tax_table
  expect_identical(tax_table(ps2), tt2)
  ## applied to tax_table object, reorder = FALSE
  tax3 <- tax_glom(tax_table(ps0), "Phylum", reorder = FALSE)
  # expected tax table
  tt3 <- tb %>%
    dplyr::arrange(first_row) %>%
    dplyr::select(first_otu, rank_names(ps0)) %>%
    tibble::column_to_rownames("first_otu") %>%
    as("matrix") %>%
    tax_table
  expect_identical(tax3, tt3)
  ## applied to tax_table object, reorder = TRUE
  tax4 <- tax_glom(tax_table(ps0), "Phylum", reorder = TRUE)
  # expected tax table
  tt4 <- tb %>%
    dplyr::arrange(Kingdom, Phylum) %>%
    dplyr::select(first_otu, rank_names(ps0)) %>%
    tibble::column_to_rownames("first_otu") %>%
    as("matrix") %>%
    tax_table
  expect_identical(tax4, tt4)
})

# TODO: test w/ non-default bad_empty being used
# TODO: test w/ tax tables as characters and as factors

# phyloseq tests
# https://github.com/joey711/phyloseq/blob/master/tests/testthat/test-merge.R
GP.chl = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
test_that("the tax_table slot is identical whether tax_glom()ed by itself or as component", {
	expect_is(tax_glom(tax_table(GP.chl), "Family"), ("taxonomyTable"))
	expect_is(n1<-tax_glom(GP.chl, "Family"), ("phyloseq"))
	expect_equal(ntaxa(n1), (4L))
	expect_equivalent(
		tax_glom(tax_table(GP.chl), taxrank="Family"),
		tax_table(tax_glom(GP.chl, taxrank="Family"))
	)
	n1 = as(tax_glom(tax_table(GP.chl), taxrank="Family", NArm=FALSE), "matrix")[, "Family"]
	n2 = tax_glom(GP.chl, taxrank="Family", NArm=FALSE)
  expect_true(setequal(n1, as(tax_table(n2), "matrix")[, "Family"]))
	expect_equal(ntaxa(n2), (5L))	
})
test_that("tax_glom() handles clearly agglomeration to one taxa", {
	expect_warning(n1 <- tax_glom(GP.chl, "Phylum"))
	expect_is(n1, ("phyloseq"))
	expect_equal(ntaxa(n1), (1L))
	expect_is(access(n1, "phy_tree"), ("NULL"))
})
test_that("tax_glom() can handle even the highest rank glom", {
  expect_warning(tax_glom(GP.chl, "Kingdom"))
  gpk = tax_glom(GlobalPatterns, "Kingdom")
  expect_is(gpk, "phyloseq")
  expect_equivalent(ntaxa(gpk), 2)
  expect_equivalent(taxa_sums(gpk), c(195598, 28021080))
})

# tip_glom--------------------------------------------------------------------- 

test_that("tip_glom() is equivalent with phyloseq::tip_glom()", {
  ps1 <- phyloseq::tip_glom(ps, 0.1)
  ps2 <- tip_glom(ps, 0.1)
  expect_identical(ps1, ps2)
})

# tree_glom-------------------------------------------------------------------- 

data(GlobalPatterns)
# ps <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
# # Prune to just the bigger clade
# phy_tree(ps) %>%
#   ggtree::ggtree +
#   ggtree::geom_nodelab()
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

# # Will test `tree_glom()` against manual merging calc'd from these branch
# # lengths:
# phy_tree(ps) %>%
#   ggtree::ggtree() +
#   ggtree::geom_tiplab() +
#   ggtree::geom_label(aes(x = branch, label = branch.length %>% round(4)))
# ggsave("/tmp/tree.pdf")

# Taxa names in order
# taxa <- c("584073", "2935", "217851", "89521", "249365", "25769", "152689",
#   "239522", "253897", "535088", "171324", "552935")

test_that("tree_glom() agrees with manual merging", {
  # Resolution = 0.001; no merging should be done
  ps1 <- tree_glom(ps, 0.001)
  expect_identical(ps1, ps)
  # Resolution = 0.03; several pairs should be merged
  ps1 <- tree_glom(ps, 0.03)
  ps2 <- ps %>%
    merge_taxa(c("253897", "239522")) %>%
    merge_taxa(c("25769", "249365")) %>%
    merge_taxa(c("89521", "217851"))
  expect_identical(ps1, ps2)
  # Resolution = 0.05; some larger groups merged
  ps1 <- tree_glom(ps, 0.05)
  ps2 <- ps %>%
    merge_taxa(c("253897", "239522", "152689")) %>%
    merge_taxa(c("25769", "249365")) %>%
    merge_taxa(c("89521", "217851", "2935")) %>%
    merge_taxa(c("552935", "171324"))
  expect_identical(ps1, ps2)
  # Try alt. merging criterion
  ps1 <- tree_glom(ps, 0.070, criterion = "max_tip_pair_dist")
  ps2 <- ps %>%
    merge_taxa(c("253897", "239522", "152689")) %>%
    merge_taxa(c("25769", "249365")) %>%
    merge_taxa(c("89521", "217851")) %>%
    merge_taxa(c("552935", "171324"))
  expect_identical(ps1, ps2)
})
