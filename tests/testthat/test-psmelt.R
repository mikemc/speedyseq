# devtools::test_file(here::here("tests", "testthat", "test-psmelt.R"))

library(dplyr)

# Test on a subset of GlobalPatterns
data(GlobalPatterns)
set.seed(20190421)
ps <- GlobalPatterns %>%
  {prune_taxa(sample(taxa_names(.), 40), .)} %>%
  tax_glom("Genus")

test_that("psmelt() functionally matches phyloseq::psmelt()", {
  # In this example, phyloseq::psmelt() drops the "Species" column from the
  # tax_table (which is full of NAs). But speedyseq::psmelt() does not
  # discard columns even if they contain no non-missing data. Other
  # differences: phyloseq's output has rownames (which seem not meaningful);
  # speedyseq's output does not; both data frames are sorted by Abundance,
  # but the row order differs in cases where Abundance and OTU names differ.
  options(stringsAsFactors = FALSE)
  tb1 <- phyloseq::psmelt(ps)
  tb2 <- psmelt(ps)
  expect_is(tb2$Kingdom, "character")
  expect_true(all_equal(tb1, select(tb2, -Species),
    ignore_col_order = FALSE, ignore_row_order = TRUE))
  expect_equal(tb1$Abundance, tb2$Abundance)
  expect_identical(tb1$OTU, tb2$OTU)
})

test_that("psmelt() provides the requested class", {
  dt <- psmelt(ps, as = "data.table")
  df <- psmelt(ps, as = "data.frame")
  tb <- psmelt(ps, as = "tibble")
  expect_is(dt, "data.table")
  expect_is(df, "data.frame")
  expect_is(tb, "tbl_df")
  expect_identical(df, dt %>% as("data.frame"))
  expect_identical(tb, df %>% as_tibble)
  # The dt -> tibble direct conversion has the ".internal.selfref" attribute
  expect_equivalent(tb, dt %>% as_tibble)
  # By default, gets `as` from "speedyseq.psmelt_class" option
  dt <- withr::with_options(
    list(speedyseq.psmelt_class = "data.table"), 
    psmelt(ps)
  )
  expect_is(dt, "data.table")
})

# Test code taken from phyloseq -----------------------------------------------
# https://github.com/joey711/phyloseq/blob/master/tests/testthat/test-plot.R

test_that("psmelt properly protects against various name collisions", {
  data("GlobalPatterns")
  gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
  ps1 = NULL
  gp1 = gp.ch
  # type-1a conflict, Abundance
  sample_data(gp1)$Abundance <- paste0("Sa-", 1:nsamples(gp1))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_Abundance" %in% colnames(ps1))
  # A different type-1a conflict, OTU
  ps1 = NULL
  gp1 = gp.ch
  sample_data(gp1)$OTU <- paste0("Sa-", 1:nsamples(gp1))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_OTU" %in% colnames(ps1))
  # A different type-1a conflict, Sample
  ps1 = NULL
  gp1 = gp.ch
  sample_data(gp1)$Sample <- paste0("Sa-", 1:nsamples(gp1))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_Sample" %in% colnames(ps1))  
  # type-1b conflict. rank_names conflict with special variables
  ps1 = NULL
  gp1 = gp.ch
  tax_table(gp1) <- cbind(tax_table(gp1), Sample=paste0("ta", taxa_names(gp1)))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("taxa_Sample" %in% colnames(ps1))  
  # type-2 conflict. Variable collision between rank_names and sample_data
  ps1 = NULL
  gp1 = gp.ch
  tax_table(gp1) <- cbind(tax_table(gp1), Primer=paste0("ta", taxa_names(gp1)))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 18L))
  expect_true("sample_Primer" %in% colnames(ps1)) 
  # All conflict types at once.
  ps1 = NULL
  gp1 = gp.ch
  sample_data(gp1)$Abundance <- paste0("Sa-", 1:nsamples(gp1))  
  sample_data(gp1)$OTU <- paste0("Sa-", 1:nsamples(gp1))
  sample_data(gp1)$Sample <- paste0("Sa-", 1:nsamples(gp1))
  tax_table(gp1) <- cbind(tax_table(gp1), Sample=paste0("ta", taxa_names(gp1)))
  tax_table(gp1) <- cbind(tax_table(gp1), Primer=paste0("ta", taxa_names(gp1)))
  expect_warning(ps1 <- psmelt(gp1))
  expect_equal(colnames(ps1)[1:3], c("OTU", "Sample", "Abundance"))
  expect_equal(dim(ps1), c(546L, 22L))
  newvars = c("sample_OTU", "sample_Sample", "sample_Abundance",
              "sample_Primer", "taxa_Sample")
  expect_true(all(newvars %in% colnames(ps1)))   
})

test_that("psmelt correctly handles phyloseq data with NULL components, and OTU tables", {
  data("GlobalPatterns")
  GP = prune_taxa(names(sort(taxa_sums(GlobalPatterns), TRUE)[1:50]),
    GlobalPatterns)
  # The objects with NULL components
  GPS = phyloseq(otu_table(GP), sample_data(GP), phy_tree(GP))
  GPT = phyloseq(otu_table(GP), tax_table(GP), phy_tree(GP))
  GPTr = phyloseq(otu_table(GP), phy_tree(GP))
  GPN = otu_table(GP)
  # Try psmelt directly. Should be no errors or warnings.
  expect_is((testT <- psmelt(GPT)), "data.frame")
  expect_is((testS <- psmelt(GPS)), "data.frame")
  expect_is((testTr <- psmelt(GPTr)), "data.frame")
  expect_is((testN <- psmelt(GPN)), "data.frame")
  # Test values of the results.
  expect_is(testT$Abundance, "numeric")
  expect_is(testT$OTU, "character")
  expect_is(testT$Sample, "character")
  expect_equivalent(
    colnames(testT),
    c("OTU", "Sample", "Abundance", "Kingdom", "Phylum", "Class", "Order",
      "Family", "Genus", "Species")
  )
  expect_equivalent(
    colnames(testS), 
    # NOTE: phyloseq instead expects c("Sample", "OTU", ... for this case
    c("OTU", "Sample",
      "Abundance", "X.SampleID", "Primer", "Final_Barcode",
      "Barcode_truncated_plus_T", "Barcode_full_length", "SampleType",
      "Description")
  )
  # Try psmelt via plot function that relies on it
  expect_is(pS <- plot_tree(GPS, color="SampleType"), "ggplot")
  expect_is(pT <- plot_tree(GPT, shape="Kingdom"), "ggplot")
  expect_is(pTr <- plot_tree(GPTr), "ggplot")
  expect_is(pN <- plot_bar(GPN), "ggplot")
  # Note, these calls produce graphical output
  expect_is((prPS <- print(pS)), "gg")
  expect_is((prPT <- print(pT)), "gg")
  expect_is((prPTr <- print(pTr)), "gg")
  expect_is((prPN <- print(pN)), "gg")
})

test_that("psmelt doesn't break when the number of taxa is 1", {
  data(GlobalPatterns)
  # tree removal warning when prune to 1 OTU.
  expect_warning(GP1 <- prune_taxa(taxa_names(GlobalPatterns)[1], GlobalPatterns))
  expect_equal(ntaxa(GP1), 1)
  df <- psmelt(GP1)
  expect_is(df, 'data.frame')
  reqnames = c("OTU", "Sample", "Abundance", "SampleType", "Kingdom", "Phylum")
  expect_true(all(reqnames %in% names(df)))
  expect_equivalent(sum(df$Abundance, na.rm = TRUE), taxa_sums(GP1))
})
