# run just this file:
# devtools::test_active_file(here::here("tests", "testthat", "test-dplyr-verbs.R"))

data(GlobalPatterns)
ps <- GlobalPatterns %>% prune_taxa(head(taxa_names(.), 10), .)

test_that("mutate and filter tax table", {
  ps1 <- ps %>%
    mutate_tax_table(
      .otu = stringr::str_c(Kingdom, .otu),
      across(Class, ~stringr::str_c("c__", .))
    ) %>%
    filter_tax_table(stringr::str_detect(Class, "^c__Th"))
  tax1 <- ps %>%
    tax_table %>%
    mutate_tax_table(
      .otu = stringr::str_c(Kingdom, .otu),
      across(Class, ~stringr::str_c("c__", .))
    ) %>%
    filter_tax_table(stringr::str_detect(Class, "^c__Th"))
  expect_identical(ntaxa(ps1), 3L)
  expect_identical(tax_table(ps1), tax1)
})

test_that("mutate and filter sample data", {
  ps1 <- ps %>%
    mutate_sample_data(
      .sample = stringr::str_c(SampleType, .sample),
      new_var = dplyr::row_number()
    ) %>%
    filter_sample_data(SampleType %in% c("Feces", "Soil"))
  sam1 <- ps %>%
    sample_data %>%
    mutate_sample_data(
      .sample = stringr::str_c(SampleType, .sample),
      new_var = dplyr::row_number()
    ) %>%
    filter_sample_data(SampleType %in% c("Feces", "Soil"))
  expect_identical(nsamples(ps1), 7L)
  expect_identical(sample_data(ps1), sam1)
})

test_that("can rename columns in taxonomy table and sample data", {
  ps1 <- GlobalPatterns %>%
    rename_tax_table(Domain = Kingdom) %>%
    rename_with_tax_table(stringr::str_to_lower) %>%
    rename_with_sample_data(janitor::make_clean_names) %>%
    rename_sample_data(sample_id = x_sample_id)

  nms <- rank_names(GlobalPatterns) %>% stringr::str_to_lower()
  nms[1] <- "domain"
  expect_identical(ps1 %>% rank_names, nms)

  nms <- sample_variables(GlobalPatterns) %>% janitor::make_clean_names()
  nms[1] <- "sample_id"
  expect_identical(ps1 %>% sample_variables, nms)
})
