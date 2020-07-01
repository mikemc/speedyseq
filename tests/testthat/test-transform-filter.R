# run just this file:
# devtools::test_file(here::here("tests", "testthat", "test-transform-filter.R"))
data(GlobalPatterns)

test_that("transform and filter functions agree with phyloseq", {
  gp1 <- filter_taxa2(GlobalPatterns, ~ sum(. > 0) > 5)
  gp2 <- phyloseq::filter_taxa(GlobalPatterns, function(x) sum(x > 0) > 5, 
    prune = TRUE)
  expect_identical(gp1, gp2)

  gp3 <- transform_sample_counts(gp1, ~ . / sum(.))
  gp4 <- phyloseq::transform_sample_counts(gp1, function(x) x / sum(x))
  expect_identical(gp3, gp4)
})
