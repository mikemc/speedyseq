# devtools::test_file(here::here("tests", "testthat", "test-merge_taxa_vec.R"))

library(magrittr)

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

test_that("merge_taxa_vec() passes basic checks", {
  # Nothing should be merged
  ps1 <- merge_taxa_vec(ps, group = taxa_names(ps))
  expect_identical(ps, ps1)
  # Should fail with an error
  expect_error(merge_taxa_vec(sample_data(ps), taxa_names(ps)))
  # Should only be one taxon; tree will be dropped
  expect_warning(ps1 <- merge_taxa_vec(ps, rep(1, ntaxa(ps))))
  expect_identical(ntaxa(ps1), 1L)
  expect_null(access(ps1, "phy_tree"))
})

test_that("merge_taxa_vec() agrees with merge_taxa", {
  taxa <- taxa_names(ps)
  group <- c(1, 9, 3, 4, 5, 3, 7, 8, 9, 9, 9, 12)
  f <- function(x) {
    x %>%
      merge_taxa(taxa[group == 9]) %>%
      merge_taxa(taxa[group == 3])
  }
  expect_identical(merge_taxa_vec(ps, group), f(ps))
  # Check correct on individual components
  x <- tax_table(ps); expect_identical(merge_taxa_vec(x, group), f(x))
  x <- phy_tree(ps);  expect_identical(merge_taxa_vec(x, group), f(x))
  x <- refseq(ps);    expect_identical(merge_taxa_vec(x, group), f(x))
  # For the otu table, allow for different taxa order
  x <- otu_table(ps)
  y1 <- merge_taxa_vec(x, group)
  y2 <- f(x)
  expect_true(setequal(taxa_names(y1), taxa_names(y2)))
  expect_identical(select_taxa(y1, taxa_names(y2)), y2)
})

test_that("merge_taxa_vec() works in both taxa orientations, with numerics and factors, and warns with NAs", {
  group <- c(1, 9, 3, 4, 5, 3, NA, 8, 9, NA, 9, 12)
  taxa <- taxa_names(ps)
  expect_warning(ps1 <- merge_taxa_vec(ps, group))
  expect_warning(ps2 <- merge_taxa_vec(t(ps), as.factor(group)) %>% t)
  ps3 <- ps %>% 
    prune_taxa(taxa[!is.na(group)], .) %>%
    merge_taxa(taxa[group %in% 3]) %>%
    merge_taxa(taxa[group %in% 9])
  expect_identical(ps1, ps2)
  expect_identical(ps1, ps3)
})

# Test results with differet tax_adjust and reorder options -------------------
#
test_that("merge_taxa_vec() gives expected results with different tax_adjust and reorder options", {
  # New test data
  idx <- c(1, 2, 5000, 5001, 10000, 18000)
  taxa <- taxa_names(GlobalPatterns)[idx]
  ps <- select_taxa(GlobalPatterns, taxa)
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
  # Add NA at intermediate rank for testing tax_adjust behavior
  tax_table(ps)[1:2, 2] <- NA
  # cbind(sum = taxa_sums(ps), tax_table(ps))
  #>        sum    Kingdom    Phylum           Class                 Order          
  #> 549322 "259"  "Archaea"  NA               "Thermoprotei"        NA             
  #> 522457 "8"    "Archaea"  NA               "Thermoprotei"        NA             
  #> 518474 "9"    "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Legionellales"
  #> 244967 "13"   "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Legionellales"
  #> 362168 "101"  "Bacteria" "Bacteroidetes"  "Bacteroidia"         "Bacteroidales"
  #> 367055 "9344" "Bacteria" "Firmicutes"     "Clostridia"          "Clostridiales"
  #>        Family            Genus          Species
  #> 549322 NA                NA             NA     
  #> 522457 NA                NA             NA     
  #> 518474 "Legionellaceae"  "Legionella"   NA     
  #> 244967 "Legionellaceae"  "Legionella"   NA     
  #> 362168 "Bacteroidaceae"  "Bacteroides"  NA     
  #> 367055 "Ruminococcaceae" "Ruminococcus" NA     
  # grouping for tests
  group <- c(3, 3, 2, 2, 1, 1)
  # First taxon in each group
  first_names <- c("549322", "518474", "362168")
  # Most abundant taxon in each group
  max_names <- c("549322", "244967", "367055")

  #### Without taxa reordering
  ## Applied to a tax table
  # The expected results for different tax_adjust options with reorder = FALSE
  # and applied to a tax_table
  tax0 <- tax_table(ps)[first_names,]
  tax1 <- tax2 <- tax0
  tax1[3, 2:7] <- tax2[3, 2:7] <- NA
  tax2[1, 3] <- NA
  expect_equal(tax0,
    merge_taxa_vec(ps %>% tax_table, group, tax_adjust = 0, reorder = FALSE))
  expect_equal(tax1,
    merge_taxa_vec(ps %>% tax_table, group, tax_adjust = 1, reorder = FALSE))
  expect_equal(tax2,
    merge_taxa_vec(ps %>% tax_table, group, tax_adjust = 2, reorder = FALSE))
  ## Applied to a phyloseq object
  # The expected results for different tax_adjust options with reorder = FALSE
  # and applied to a phyloseq object
  tax0 <- tax_table(ps)[max_names,]
  tax1 <- tax2 <- tax0
  tax1[3, 2:7] <- tax2[3, 2:7] <- NA
  tax2[1, 3] <- NA
  expected_otu_mat <- rowsum(otu_table(ps), group, reorder = FALSE)
  rownames(expected_otu_mat) <- max_names
  ps0 <- merge_taxa_vec(ps, group, tax_adjust = 0, reorder = FALSE)
  ps1 <- merge_taxa_vec(ps, group, tax_adjust = 1, reorder = FALSE)
  ps2 <- merge_taxa_vec(ps, group, tax_adjust = 2, reorder = FALSE)
  expect_identical(tax0, ps0 %>% tax_table)
  expect_identical(tax1, ps1 %>% tax_table)
  expect_identical(tax2, ps2 %>% tax_table)
  expect_identical(expected_otu_mat, ps0 %>% otu_table %>% as("matrix"))
  expect_identical(expected_otu_mat, ps1 %>% otu_table %>% as("matrix"))
  expect_identical(expected_otu_mat, ps2 %>% otu_table %>% as("matrix"))
  expect_identical(refseq(ps)[max_names], ps0 %>% refseq)
  expect_identical(refseq(ps)[max_names], ps1 %>% refseq)
  expect_identical(refseq(ps)[max_names], ps2 %>% refseq)
  expect_identical(ape::keep.tip(phy_tree(ps), max_names), ps0 %>% phy_tree)
  expect_identical(ape::keep.tip(phy_tree(ps), max_names), ps1 %>% phy_tree)
  expect_identical(ape::keep.tip(phy_tree(ps), max_names), ps2 %>% phy_tree)

  #### With taxa reordering
  ps@phy_tree <- NULL
  ## Applied to a phyloseq object
  # The expected results for different tax_adjust options with reorder = TRUE
  # and applied to a phyloseq object with no tree
  tax0 <- tax_table(ps)[rev(max_names),]
  tax1 <- tax2 <- tax0
  tax1[1, 2:7] <- tax2[1, 2:7] <- NA
  tax2[3, 3] <- NA
  expected_otu_mat <- rowsum(otu_table(ps), group, reorder = TRUE)
  rownames(expected_otu_mat) <- rev(max_names)
  ps0 <- merge_taxa_vec(ps, group, tax_adjust = 0, reorder = TRUE)
  ps1 <- merge_taxa_vec(ps, group, tax_adjust = 1, reorder = TRUE)
  ps2 <- merge_taxa_vec(ps, group, tax_adjust = 2, reorder = TRUE)
  expect_identical(tax0, ps0 %>% tax_table)
  expect_identical(tax1, ps1 %>% tax_table)
  expect_identical(tax2, ps2 %>% tax_table)
  expect_identical(expected_otu_mat, ps0 %>% otu_table %>% as("matrix"))
  expect_identical(expected_otu_mat, ps1 %>% otu_table %>% as("matrix"))
  expect_identical(expected_otu_mat, ps2 %>% otu_table %>% as("matrix"))
  expect_identical(refseq(ps)[rev(max_names)], ps0 %>% refseq)
  expect_identical(refseq(ps)[rev(max_names)], ps1 %>% refseq)
  expect_identical(refseq(ps)[rev(max_names)], ps2 %>% refseq)
  ## Applied to a tax table
  # The expected results for different tax_adjust options with reorder = TRUE
  # and applied to a tax_table
  tax0 <- tax_table(ps)[rev(first_names),]
  tax1 <- tax2 <- tax0
  tax1[1, 2:7] <- tax2[1, 2:7] <- NA
  tax2[3, 3] <- NA
  expect_equal(tax0,
    merge_taxa_vec(ps %>% tax_table, group, tax_adjust = 0, reorder = TRUE))
  expect_equal(tax1,
    merge_taxa_vec(ps %>% tax_table, group, tax_adjust = 1, reorder = TRUE))
  expect_equal(tax2,
    merge_taxa_vec(ps %>% tax_table, group, tax_adjust = 2, reorder = TRUE))

  #### applied to just refseq objects
  expect_equal(
    refseq(ps)[first_names],
    merge_taxa_vec(ps %>% refseq, group, reorder = FALSE)
  )
  expect_equal(
    refseq(ps)[rev(first_names)],
    merge_taxa_vec(ps %>% refseq, group, reorder = TRUE)
  )
})
