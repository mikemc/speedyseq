context("Equivalence of tax_glom() with phyloseq::tax_glom()")

library(dplyr)

data(GlobalPatterns)

# Pick a subset to speed up the tests
set.seed(20190421)
ps <- GlobalPatterns %>%
    {prune_taxa(sample(taxa_names(.), 1e3), .)}

# tax_glom--------------------------------------------------------------------- 

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for an intermediate rank with default arguments", {
    ps1 <- phyloseq::tax_glom(ps, "Phylum")
    ps2 <- tax_glom(ps, "Phylum")
    expect_equal(taxa_names(ps1), taxa_names(ps2))
    expect_equal(otu_table(ps1), otu_table(ps2))
    expect_equal(sample_data(ps1), sample_data(ps2))
    expect_equal(tax_table(ps1), tax_table(ps2))
    expect_equal(phy_tree(ps1), phy_tree(ps2))
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for an intermediate rank with NArm = FALSE", {
    ps1 <- phyloseq::tax_glom(ps, "Family", NArm = FALSE)
    ps2 <- tax_glom(ps, "Family", NArm = FALSE)
    expect_equal(taxa_names(ps1), taxa_names(ps2))
    expect_equal(otu_table(ps1), otu_table(ps2))
    expect_equal(sample_data(ps1), sample_data(ps2))
    expect_equal(tax_table(ps1), tax_table(ps2))
    expect_equal(phy_tree(ps1), phy_tree(ps2))
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for the highest rank", {
    ps1 <- phyloseq::tax_glom(ps, "Kingdom")
    ps2 <- tax_glom(ps, "Kingdom")
    expect_equal(taxa_names(ps1), taxa_names(ps2))
    expect_equal(otu_table(ps1), otu_table(ps2))
    expect_equal(sample_data(ps1), sample_data(ps2))
    expect_equal(tax_table(ps1), tax_table(ps2))
    expect_equal(phy_tree(ps1), phy_tree(ps2))
})

test_that("tax_glom() is equivalent with phyloseq::tax_glom() for the lowest rank", {
    ps1 <- phyloseq::tax_glom(ps, "Species")
    ps2 <- tax_glom(ps, "Species")
    expect_equal(taxa_names(ps1), taxa_names(ps2))
    expect_equal(otu_table(ps1), otu_table(ps2))
    expect_equal(sample_data(ps1), sample_data(ps2))
    expect_equal(tax_table(ps1), tax_table(ps2))
    expect_equal(phy_tree(ps1), phy_tree(ps2))
})


# TODO: test w/ non-default bad_empty being used
# TODO: test w/ tax tables as characters and as factors
# j
#

# tip_glom--------------------------------------------------------------------- 

# Pick a subset to speed up the tests
set.seed(20190421)
ps <- GlobalPatterns %>%
    {prune_taxa(sample(taxa_names(.), 2e2), .)}

test_that("tip_glom() is equivalent with phyloseq::tip_glom()", {
    ps1 <- phyloseq::tip_glom(ps, 0.1)
    ps2 <- tip_glom(ps, 0.1)
    expect_equal(taxa_names(ps1), taxa_names(ps2))
    expect_equal(otu_table(ps1), otu_table(ps2))
    expect_equal(sample_data(ps1), sample_data(ps2))
    expect_equal(tax_table(ps1), tax_table(ps2))
    expect_equal(phy_tree(ps1), phy_tree(ps2))
})
