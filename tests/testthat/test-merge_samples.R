# run just this file:
# devtools::test_active_file(here::here("tests", "testthat", "test-merge_samples.R"))

# helpers ---------------------------------------------------------------------

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

test_that("`sample_data_stable()` maintains names", {
  x <- data.frame(var1 = letters[1:3], var2 = 7:9)
  expect_identical(sample_data_stable(x) %>% sample_names, as.character(1:3))
  rownames(x) <- c("3", "two", "1")
  expect_identical(sample_data_stable(x) %>% sample_names, c("3", "two", "1"))
})

# merge_samples2 --------------------------------------------------------------

test_that("Test `merge_samples2()` on the `enterotype` dataset", {
  data(enterotype)
  ps <- enterotype %>% prune_taxa(head(taxa_names(.), 10), .)
  sample_data(ps) <- sample_data(ps) %>%
    transform(Project.ClinicalStatus = Project:ClinicalStatus)
  # Will get a warning since `group` has NAs
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
  # Should also work on columns that are numbers, or factors whose levels are
  # numbers, even if the values are of the form 1:n. (Issue #52)
  expect_warning(ps2 <- merge_samples2(ps, "Enterotype"))
  sample_data(ps)$Enterotype <- sample_data(ps)$Enterotype %>% as.integer
  expect_warning(ps3 <- merge_samples2(ps, "Enterotype"))
  sample_data(ps3)$Enterotype <- sample_data(ps3)$Enterotype %>% as.factor
  expect_identical(ps2, ps3)
  # Test with alternate abundance-merging function
  expect_warning(
    ps0 <- merge_samples2(ps, "Project.ClinicalStatus",
      fun_otu = mean,
      funs = list(Age = mean)
    )
  )
})

test_that("Test that group names are preserved", {
  # For Issue #52
  sam <- tibble::tibble(sample_id = letters[1:3], group_var = 1:3) %>%
    sample_data
  x <- merge_samples2(sam, "group_var")
  expect_identical(sample_names(x), as.character(sam$group_var))
})
