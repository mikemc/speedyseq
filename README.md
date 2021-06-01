
<!-- README.md is generated from README.Rmd. Please edit that file -->

# speedyseq

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/179732395.svg)](https://zenodo.org/badge/latestdoi/179732395)
[![Travis build
status](https://travis-ci.com/mikemc/speedyseq.svg?branch=main)](https://travis-ci.com/mikemc/speedyseq)
[![Codecov test
coverage](https://codecov.io/gh/mikemc/speedyseq/branch/main/graph/badge.svg)](https://codecov.io/gh/mikemc/speedyseq?branch=main)
<!-- badges: end -->

Speedyseq is an R package for microbiome data analysis that extends the
popular [phyloseq](https://joey711.github.io/phyloseq/) package.
Speedyseq began with the limited goal of providing faster versions of
phyloseq’s plotting and taxonomic merging functions, but now contains a
growing number of enhancements to phyloseq which I have found useful.

## Installation

Install the current development version with the remotes package,

``` r
# install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
```

## Usage

Method 1: Call speedyseq functions explicitly when you want to use
speedyseq’s version instead of phyloseq. This method ensures that you do
not unintentionally call speedyseq’s version of a phyloseq function.

``` r
library(phyloseq)
data(GlobalPatterns)
system.time(
  # Calls phyloseq's psmelt
  df1 <- psmelt(GlobalPatterns) # slow
)
#>    user  system elapsed 
#>   6.320   0.063   6.390
system.time(
  df2 <- speedyseq::psmelt(GlobalPatterns) # fast
)
#>    user  system elapsed 
#>   0.344   0.004   0.245
dplyr::all_equal(df1, df2, ignore_row_order = TRUE)
#> [1] TRUE
detach(package:phyloseq)
```

Method 2: Load speedyseq, which will load phyloseq and all speedyseq
functions and cause calls to the overlapping function names to go to
speedyseq by default.

``` r
library(speedyseq)
#> Loading required package: phyloseq
#> 
#> Attaching package: 'speedyseq'
#> The following objects are masked from 'package:phyloseq':
#> 
#>     filter_taxa, plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom, tip_glom,
#>     transform_sample_counts
data(GlobalPatterns)
system.time(
  ps1 <- phyloseq::tax_glom(GlobalPatterns, "Genus") # slow
)
#>    user  system elapsed 
#>  31.856   0.106  32.031
system.time(
  # Calls speedyseq's tax_glom
  ps2 <- tax_glom(GlobalPatterns, "Genus") # fast
)
#>    user  system elapsed 
#>   0.241   0.000   0.230
```

Loading speedyseq will also load the
[magrittr](https://magrittr.tidyverse.org/) pipe (`%>%`) to allow pipe
chains with phyloseq objects,

``` r
gp.filt.prop <- GlobalPatterns %>%
  filter_taxa2(~ sum(. > 0) > 5) %>%
  transform_sample_counts(~ . / sum(.))
```

## Features

### Faster implementations of phyloseq functions

-   `psmelt()` and the plotting functions that use it: `plot_bar()`,
    `plot_heatmap()`, and `plot_tree()`.
-   The taxonomic merging functions `tax_glom()` and `tip_glom()`.
    Speedyseq’s `tip_glom()` also has significantly lower memory usage.

These functions should generally function as drop-in replacements for
phyloseq’s versions, with additional arguments allowing for modified
behavior. Differences in row order (for `psmelt()`) and taxon order (for
`tax_glom()`) can occur; see
[Changelog](https://mikemc.github.io/speedyseq/news/index.html) for
details.

### New taxonomic merging functions

-   A general-purpose merging function `merge_taxa_vec()` that provides
    a vectorized version of phyloseq’s `merge_taxa()` function.
-   A function `tree_glom()` that performs direct phylogenetic merging
    of taxa. This function provides an alternative to the indirect
    phylogenetic merging done by `tip_glom()` that is much faster and
    arguably more intuitive.

See the [Changelog](https://mikemc.github.io/speedyseq/news/index.html)
for details and examples.

### Enhancements and additions to other phyloseq functions

See the [online
documentation](https://mikemc.github.io/speedyseq/reference/index.html)
for an up-to-date list and usage information and the
[Changelog](https://mikemc.github.io/speedyseq/news/index.html) for
further information.
