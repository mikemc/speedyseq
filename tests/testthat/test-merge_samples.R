# run just this file:
# devtools::test_file(here::here("tests", "testthat", "test-merge_samples.R"))

# unique_or_na ----------------------------------------------------------------

test_that("unique_or_na() returns expected results", {
  ## Atomic vectors
  # Characters
  expect_identical(unique_or_na(NA_character_), NA_character_)
  expect_identical(unique_or_na("a"), "a")
  expect_identical(unique_or_na(c("a", "b")), NA_character_)
  expect_identical(unique_or_na(c("a", NA_character_)), NA_character_)
  # Integers
  expect_identical(unique_or_na(NA_integer_), NA_integer_)
  expect_identical(unique_or_na(1L), 1L)
  expect_identical(unique_or_na(c(1L, 2L)), NA_integer_)
  expect_identical(unique_or_na(c(1L, NA_integer_)), NA_integer_)
  # Doubles
  expect_identical(unique_or_na(NA_real_), NA_real_)
  expect_identical(unique_or_na(1), 1)
  expect_identical(unique_or_na(c(1, 2)), NA_real_)
  expect_identical(unique_or_na(c(1, NA_real_)), NA_real_)
  # Booleans
  expect_identical(unique_or_na(NA), NA)
  expect_identical(unique_or_na(TRUE), TRUE)
  expect_identical(unique_or_na(c(TRUE, FALSE)), NA)
  expect_identical(unique_or_na(c(TRUE, NA)), NA)
  # Factors
  expect_identical(
    unique_or_na(factor(letters[1:4], levels = letters[4:1])),
    factor(NA, levels = letters[4:1])
  )
  expect_identical(
    unique_or_na(factor(c("a", "a"), levels = letters[4:1])),
    factor("a", levels = letters[4:1])
  )
  expect_identical(
    unique_or_na(factor(c("a", "a"), levels = letters[4:1], 
      labels = paste("letter", letters[4:1]))),
    factor("a", levels = letters[4:1], labels = paste("letter", letters[4:1]))
  )
  expect_identical(
    unique_or_na(factor(letters[1:4], levels = letters[4:1], ordered = TRUE)),
    factor(NA, levels = letters[4:1], ordered = TRUE)
  )
  ## Lists
  expect_identical(unique_or_na(list(1,2)), NA)
  expect_identical(unique_or_na(list(1,1)), 1)
})

# merge_samples2 --------------------------------------------------------------

test_that("Test `merge_samples2()` on `enterotype` dataset", {
  data(enterotype)
  ps <- enterotype 
  sample_data(ps) <- sample_data(ps) %>%
    transform(Project.ClinicalStatus = Project:ClinicalStatus)
  expect_warning(
    ps0 <- merge_samples2(ps, "Project.ClinicalStatus", 
      funs = list(Age = mean))
  )
  # Types of the sample vars should be unchanged
  purrr::walk2(sample_data(ps), sample_data(ps0), 
    ~expect_equal(typeof(.x), typeof(.y))
  )
  # Check OTU table, including that order is kept
  grp <- sample_data(ps)$Project.ClinicalStatus
  expect_warning(
    x <- rowsum(t(otu_table(ps)), grp, reorder = FALSE) %>% 
      {.[!is.na(rownames(.)),]} %>% t
  )
  expect_identical(otu_table(ps0) %>% as("matrix"), x)
  # Result should be identical if computed with taxa as cols
  expect_warning(
    ps1 <- merge_samples2(ps %>% t, "Project.ClinicalStatus", 
      funs = list(Age = mean)) %>% t
  )
  expect_identical(ps0, ps1)
  # Remaining components should be the same as original
  expect_identical(tax_table(ps0), tax_table(ps))
})
