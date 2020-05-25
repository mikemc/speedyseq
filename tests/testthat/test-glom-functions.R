context("Equivalence of tax_glom() with phyloseq::tax_glom()")

library(magrittr) # Imported for the pipe

# Function to test for equality of two phyloseq objects
compare_phyloseq <- function(x, y) {
  expect_equal( taxa_names(x),    taxa_names(y) )
  expect_equal( sample_names(x),  sample_names(y) )
  expect_equal( taxa_are_rows(x), taxa_are_rows(y) )
  expect_equal( slot(x, "otu_table"), slot(y, "otu_table") )
  expect_equal( slot(x, "sam_data"),  slot(y, "sam_data")  )
  expect_equal( slot(x, "tax_table"), slot(y, "tax_table") )
  expect_equal( slot(x, "phy_tree"),  slot(y, "phy_tree")  )
  expect_equal( slot(x, "refseq"),    slot(y, "refseq")    )
}

data(GlobalPatterns)

# Pick a subset to speed up the tests
set.seed(20190421)
ps <- GlobalPatterns %>%
    {prune_taxa(sample(taxa_names(.), 2e2), .)}
tax <- tax_table(ps)


# tax_glom--------------------------------------------------------------------- 

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for an intermediate rank with default arguments", {
    ps1 <- phyloseq::tax_glom(ps, "Phylum")
    ps2 <- tax_glom(ps, "Phylum")
    compare_phyloseq(ps1, ps2)
    # Also test on a tax_table only
    tax1 <- phyloseq::tax_glom(tax, "Phylum")
    tax2 <- tax_glom(tax, "Phylum")
    expect_equal(tax1, tax2)
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for an intermediate rank with NArm = FALSE", {
    ps1 <- phyloseq::tax_glom(ps, "Family", NArm = FALSE)
    ps2 <- tax_glom(ps, "Family", NArm = FALSE)
    compare_phyloseq(ps1, ps2)
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for the highest rank", {
    ps1 <- phyloseq::tax_glom(ps, "Kingdom")
    ps2 <- tax_glom(ps, "Kingdom")
    compare_phyloseq(ps1, ps2)
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for the lowest rank", {
    ps1 <- phyloseq::tax_glom(ps, "Species")
    ps2 <- tax_glom(ps, "Species")
    compare_phyloseq(ps1, ps2)
})


# TODO: test w/ non-default bad_empty being used
# TODO: test w/ tax tables as characters and as factors
# j
#

# tip_glom--------------------------------------------------------------------- 

test_that("tip_glom() is equivalent with phyloseq::tip_glom()", {
    ps1 <- phyloseq::tip_glom(ps, 0.1)
    ps2 <- tip_glom(ps, 0.1)
    compare_phyloseq(ps1, ps2)
})

# tree_glom-------------------------------------------------------------------- 

data(GlobalPatterns)
ps <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
# # Prune to just the bigger clade
# phy_tree(ps) %>%
#   ggtree::ggtree +
#   ggtree::geom_nodelab()
ps0 <- prune_taxa(
  ape::extract.clade(phy_tree(ps), "0.589.38") %>% .$tip.label,
  ps
)
# add reference sequences
rs <- phy_tree(ps0) %>% 
  phangorn::simSeq(10) %>%
  purrr::map(factor, levels = as.character(1:4)) %>%
  purrr::map(forcats::fct_recode, A = "1", C = "2", G = "3", T = "4") %>%
  purrr::map(as.character) %>%
  purrr::map_chr(paste, collapse = "") %>%
  Biostrings::DNAStringSet()
ps0 <- merge_phyloseq(ps0, rs)
rm(rs)

# # Will test `tree_glom()` against manual merging calc'd from these branch
# # lengths:
# phy_tree(ps0) %>%
#   ggtree::ggtree() +
#   ggtree::geom_tiplab() +
#   ggtree::geom_label(aes(x = branch, label = branch.length %>% round(4)))
# ggsave("/tmp/tree.pdf")

# Taxa names in order
taxa <- c("584073", "2935", "217851", "89521", "249365", "25769", "152689",
  "239522", "253897", "535088", "171324", "552935")

test_that("tree_glom() agrees with manual merging", {
  # Resolution = 0.001; no merging should be done
  ps1 <- tree_glom(ps0, 0.001)
  compare_phyloseq(ps1, ps0)
  # Resolution = 0.03; several pairs should be merged
  ps1 <- tree_glom(ps0, 0.03)
  ps2 <- ps0 %>%
    merge_taxa(c("253897", "239522")) %>%
    merge_taxa(c("25769", "249365")) %>%
    merge_taxa(c("89521", "217851"))
  compare_phyloseq(ps1, ps2)
  # Resolution = 0.05; some larger groups merged
  ps1 <- tree_glom(ps0, 0.05)
  ps2 <- ps0 %>%
    merge_taxa(c("253897", "239522", "152689")) %>%
    merge_taxa(c("25769", "249365")) %>%
    merge_taxa(c("89521", "217851", "2935")) %>%
    merge_taxa(c("552935", "171324"))
  compare_phyloseq(ps1, ps2)
  # Try alt. merging criterion
  ps1 <- tree_glom(ps0, 0.070, criterion = "max_tip_pair_dist")
  ps2 <- ps0 %>%
    merge_taxa(c("253897", "239522", "152689")) %>%
    merge_taxa(c("25769", "249365")) %>%
    merge_taxa(c("89521", "217851")) %>%
    merge_taxa(c("552935", "171324"))
  compare_phyloseq(ps1, ps2)
})

# merge_taxa_vec --------------------------------------------------------------

test_that("merge_taxa_vec() passes basic checks", {
  # Nothing should be merged
  ps1 <- merge_taxa_vec(ps0, group = taxa_names(ps0))
  compare_phyloseq(ps0, ps1)
  # Should fail with an error
  expect_error(merge_taxa_vec(sample_data(ps0), taxa_names(ps0)))
  # Should only be one taxon; tree will be dropped
  expect_warning(ps1 <- merge_taxa_vec(ps0, rep(1, ntaxa(ps0))))
  expect_equal(ntaxa(ps1), 1)
  expect_null(access(ps1, "phy_tree"))
})

test_that("merge_taxa_vec() agrees with merge_taxa", {
  group <- c(1, 9, 3, 4, 5, 3, 7, 8, 9, 9, 9, 12)
  f <- function(x) {
    x %>%
      merge_taxa(c("2935", "253897", "535088", "171324")) %>%
      merge_taxa(c("217851", "25769"))
  }
  compare_phyloseq(merge_taxa_vec(ps0, group), f(ps0))
  # Check correct on individual components
  x <- tax_table(ps0); expect_equal(merge_taxa_vec(x, group), f(x))
  x <- phy_tree(ps0);  expect_equal(merge_taxa_vec(x, group), f(x))
  x <- refseq(ps0);    expect_equal(merge_taxa_vec(x, group), f(x))
  # For the otu table, allow for different taxa order
  x <- otu_table(ps0)
  expect_true(
    dplyr::all_equal(
      merge_taxa_vec(x, group)@.Data %>% tibble::as_tibble(rownames = "rn"),
      f(x) %>% .@.Data %>% tibble::as_tibble(rownames = "rn")
    )
  )
})

test_that("merge_taxa_vec() works in both taxa orientations, with numerics and factors, and warns with NAs", {
  group <- c(1, 9, 3, 4, 5, 3, NA, 8, 9, NA, 9, 12)
  expect_warning(ps1 <- merge_taxa_vec(ps0, group))
  expect_warning(ps2 <- merge_taxa_vec(t(ps0), as.factor(group)) %>% t)
  ps3 <- ps0 %>% 
    prune_taxa(taxa[!is.na(group)], .) %>%
    merge_taxa(c("2935", "253897", "171324")) %>%
    merge_taxa(c("217851", "25769"))
  compare_phyloseq(ps1, ps2)
  compare_phyloseq(ps1, ps3)
})
