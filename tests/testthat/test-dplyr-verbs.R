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

test_that("can select and relocate columns in taxonomy table and sample data", {
  ps1 <- GlobalPatterns %>%
    select_tax_table(Phylum, Genus:Species) %>%
    select_sample_data(-dplyr::contains("Barcode"))
  expect_identical(c("Phylum", "Genus", "Species"), rank_names(ps1))
  nms <- GlobalPatterns %>% sample_variables %>% 
    stringr::str_subset("Barcode", negate = TRUE)
  expect_identical(nms, sample_variables(ps1))
  ps2 <- ps1 %>%
    relocate_sample_data(SampleType) %>%
    relocate_tax_table(Phylum, .after = dplyr::last_col())
  expect_identical(sample_variables(ps2),
    c("SampleType", setdiff(sample_variables(ps1), "SampleType"))
  )
  expect_identical(rank_names(ps2),
    c(setdiff(rank_names(ps1), "Phylum"), "Phylum")
  )
})

test_that("can join on sample data", {
  data(GlobalPatterns)

  ps1 <- GlobalPatterns %>%
    select_sample_data(!contains("Barcode"))
  y <- GlobalPatterns %>%
    sample_data %>%
    select_sample_data(contains("Barcode")) %>%
    as_tibble
  ps2 <- ps1 %>% left_join_sample_data(y, by = ".sample")

  z <- GlobalPatterns %>% 
    relocate_sample_data(contains("Barcode"), .after = dplyr::last_col()) %>% 
    sample_data

  expect_identical(ps2 %>% sample_data, z)
})
