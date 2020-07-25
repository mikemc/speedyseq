# run just this file:
# devtools::test_file(here::here("tests", "testthat", "test-constructors.R"))

data(GlobalPatterns)

test_that("results from tibble constructors match starting objects", {
  tax <- tax_table(GlobalPatterns) %>% 
    as("matrix") %>% 
    tibble::as_tibble(rownames = "taxon") %>%
    tax_table
  expect_identical(tax, tax_table(GlobalPatterns))

  sam <- sample_data(GlobalPatterns) %>% 
    as("data.frame") %>% 
    tibble::as_tibble(rownames = "sample") %>%
    sample_data
  expect_identical(sam, sample_data(GlobalPatterns))

  otu1 <- otu_table(GlobalPatterns) %>%
    as("matrix") %>% 
    tibble::as_tibble(rownames = "taxon") %>%
    otu_table(taxa_are_rows = TRUE)
  expect_identical(otu1, otu_table(GlobalPatterns))

  otu2 <- otu_table(GlobalPatterns) %>%
    t %>%
    as("matrix") %>% 
    tibble::as_tibble(rownames = "taxon") %>%
    otu_table(taxa_are_rows = FALSE)
  expect_identical(otu2, otu_table(GlobalPatterns) %>% t)
})
