context("Equivalence of psmelt() with phyloseq::psmelt()")

library(dplyr)

# Test on a subset of GlobalPatterns
data(GlobalPatterns)
set.seed(20190421)
ps <- GlobalPatterns %>%
    {prune_taxa(sample(taxa_names(.), 200), .)} %>%
    tax_glom("Genus")

test_that("psmelt() functionally matches phyloseq::psmelt()", {
    # In this example, phyloseq::psmelt() drops the "Species" column from the
    # tax_table (which is full of NAs). But speedyseq::psmelt() does not
    # discard columns even if they contain no non-missing data. Other
    # differences: phyloseq's output has rownames (which seem not meaningful);
    # speedyseq's output does not; both data frames are sorted by Abundance,
    # but the row order differs in cases where Abundance and OTU names differ.
    options(stringsAsFactors = TRUE)
    tb1 <- phyloseq::psmelt(ps)
    tb2 <- psmelt(ps)
    expect_true(is.factor(tb1$Kingdom))
    expect_true(all_equal(tb1, select(tb2, -Species),
        ignore_col_order = FALSE, ignore_row_order = TRUE))
    expect_equal(tb1$Abundance, tb2$Abundance)
    expect_equal(tb1$OTU, tb2$OTU)
    options(stringsAsFactors = FALSE)
    tb1 <- phyloseq::psmelt(ps)
    tb2 <- psmelt(ps)
    expect_false(is.factor(tb1$Kingdom))
    expect_true(all_equal(tb1, select(tb2, -Species),
        ignore_col_order = FALSE, ignore_row_order = TRUE))
    expect_equal(tb1$Abundance, tb2$Abundance)
    expect_equal(tb1$OTU, tb2$OTU)
})
