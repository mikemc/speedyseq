
<!-- README.md is generated from README.Rmd. Please edit that file -->

# speedyseq

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/179732395.svg)](https://zenodo.org/badge/latestdoi/179732395)
[![Travis build
status](https://travis-ci.org/mikemc/speedyseq.svg?branch=master)](https://travis-ci.org/mikemc/speedyseq)
[![Codecov test
coverage](https://codecov.io/gh/mikemc/speedyseq/branch/master/graph/badge.svg)](https://codecov.io/gh/mikemc/speedyseq?branch=master)
<!-- badges: end -->

Speedyseq aims to accelerate
[phyloseq](https://joey711.github.io/phyloseq/) operations that can be
very slow on large datasets. Current reimplementations of phyloseq
functions include

  - A faster version of phyloseq’s `psmelt()` and the plotting functions
    that make use of it (`plot_bar()`, `plot_heatmap()`, and
    `plot_tree()`).
  - Faster versions of phyloseq’s taxonomic merging functions
    `tax_glom()` and `tip_glom()`. Speedyseq’s `tip_glom()` also has
    significantly lower memory usage.

My general aim is for these functions to be drop-in replacements for
phyloseq’s versions; however, there are small differences that should
not affect most use cases. In some functions, I have added optional
arguments to allow modifying the phyloseq behavior. See
[NEWS.md](./NEWS.md) for information about these differences and
enhancements.

New functions that provide additional types of taxonomic merging include

  - A general-purpose merging function `merge_taxa_vec()` that provides
    a vectorized version of phyloseq’s `merge_taxa()` function.
  - A function `tree_glom()` that performs direct phylogenetic merging
    of taxa. This function provides an alternative to the indirect
    phylogenetic merging done by `tip_glom()` that is much faster and
    arguably more intuitive.

See [NEWS.md](./NEWS.md) for details and examples.

## Installation

Install the current development version with the remotes package,

``` r
# install.packages("remotes")
remotes::install_github("mikemc/speedyseq")
```

## Usage

Method 1: Call speedyseq functions explicitly when you want to use
speedyseq’s version instead of phyloseq:

``` r
library(phyloseq)
data(GlobalPatterns)
system.time(
  # Calls phyloseq's psmelt
  df1 <- psmelt(GlobalPatterns) # slow
)
#>    user  system elapsed 
#>  93.662   0.127  94.037
system.time(
  df2 <- speedyseq::psmelt(GlobalPatterns) # fast
)
#>    user  system elapsed 
#>   0.299   0.000   0.177
dplyr::all_equal(df1, df2, ignore_row_order = TRUE)
#> [1] TRUE
detach(package:phyloseq)
```

Method 2: Load speedyseq, which will load phyloseq and cause calls to
the overlapping function names to go to speedyseq by default:

``` r
library(speedyseq)
#> Loading required package: phyloseq
#> 
#> Attaching package: 'speedyseq'
#> The following objects are masked from 'package:phyloseq':
#> 
#>     plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom, tip_glom
data(GlobalPatterns)
system.time(
  ps1 <- phyloseq::tax_glom(GlobalPatterns, "Genus") # slow
)
#>    user  system elapsed 
#>  35.961   0.140  36.231
system.time(
  # Calls speedyseq's tax_glom
  ps2 <- tax_glom(GlobalPatterns, "Genus") # fast
)
#>    user  system elapsed 
#>   0.232   0.006   0.240
all.equal(otu_table(ps1), otu_table(ps2))
#> [1] TRUE
all.equal(tax_table(ps1), tax_table(ps2))
#> [1] TRUE
```
