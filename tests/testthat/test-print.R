data(GlobalPatterns)

test_that("print(x) returns x without error", {
  # Note, ape's print.phylo method does not return anything
  slts <- setdiff(slotNames(GlobalPatterns), "phy_tree")
  for (s in slts) {
    x <- slot(GlobalPatterns, s)
    y <- slot(GlobalPatterns, s) %>% print
    expect_identical(x, y)
  }
  # Test when only 1 column to ensure no "drop" bugs
  x <- tibble::tibble(sample = letters[1:3], var = 1:3) %>% sample_data
  expect_identical(x, print(x))
})
