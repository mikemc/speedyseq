library(dplyr)

data(GlobalPatterns)

set.seed(20190421)
ps <- GlobalPatterns %>%
  {prune_taxa(sample(taxa_names(.), 1e2), .)} %>%
  {prune_samples(sample(sample_names(.), 5), .)}

taxrank <- "Phylum"
tax_table(ps)[20:25, taxrank] <- NA
rnks <- seq(1, which(rank_names(ps) == taxrank))
group <- tax_table(ps)[, rnks] %>%
  as("matrix") %>%
  apply(1, paste, collapse = ";")
group[20:25] <- NA

test_that("merge_taxa_vec() agrees with tax_glom() with `NArm = TRUE`", {
  ps1 <- tax_glom(ps, taxrank)
  expect_warning(ps2 <- merge_taxa_vec(ps, group))
  expect_equal(taxa_names(ps1), taxa_names(ps2))
  expect_equal(otu_table(ps1), otu_table(ps2))
  expect_equal(sample_data(ps1), sample_data(ps2))
  # Tax tables should only agree to `taxrank`
  expect_equal(tax_table(ps1)[, rnks], tax_table(ps2)[, rnks])
  expect_equal(phy_tree(ps1), phy_tree(ps2))
})
