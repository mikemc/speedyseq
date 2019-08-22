<!-- README.md is generated from README.Rmd. Please edit that file -->


# speedyseq

<!-- badges: start -->
[![Travis build
status](https://travis-ci.org/mikemc/speedyseq.svg?branch=master)](https://travis-ci.org/mikemc/speedyseq)
[![Codecov test
coverage](https://codecov.io/gh/mikemc/speedyseq/branch/master/graph/badge.svg)](https://codecov.io/gh/mikemc/speedyseq?branch=master)
<!-- badges: end -->

The goal of speedyseq is to accelerate common operations that are currently
very slow in phyloseq: the function `psmelt` (and the plotting functions that
use it) and the taxonomic aggregation functions `tax_glom()` and (hopefully
eventually) `tip_glom()`.

The current version of speedyseq reimplements `psmelt()` to be faster than
phyloseq's version, and includes copies of `plot_bar()`, `plot_heatmap()`, and
`plot_tree()` so that when called from speedyseq the faster `psmelt()` will be
used. It also implements a faster version of `tax_glom()`.

## Installation

Install with devtools
```r
# install.packages("devtools")
devtools::install_github("mikemc/speedyseq")
```

## Usage

Method 1: Call speedyseq functions explicitly when you want to use
speedyseq's version instead of phyloseq:

```r
library(phyloseq)
data(GlobalPatterns)
system.time(
    # Calls phyloseq's psmelt
    df1 <- psmelt(GlobalPatterns) # slow
)
#>    user  system elapsed 
#> 101.174   0.236 102.101
system.time(
    df2 <- speedyseq::psmelt(GlobalPatterns) # fast
)
#>    user  system elapsed 
#>   0.736   0.004   0.741
dplyr::all_equal(df1, df2, ignore_row_order = TRUE)
#> [1] TRUE
detach(package:phyloseq)
```

Method 2: Load speedyseq, which will load phyloseq and cause calls to the
overlapping function names to go to speedyseq by default:

```r
library(speedyseq)
#> Loading required package: phyloseq
#> 
#> Attaching package: 'speedyseq'
#> The following objects are masked from 'package:phyloseq':
#> 
#>     plot_bar, plot_heatmap, plot_tree, psmelt, tax_glom
data(GlobalPatterns)
system.time(
    ps1 <- phyloseq::tax_glom(GlobalPatterns, "Genus") # slow
)
#>    user  system elapsed 
#>  32.827   0.089  33.013
system.time(
    # Calls speedyseq's tax_glom
    ps2 <- tax_glom(GlobalPatterns, "Genus") # fast
)
#>    user  system elapsed 
#>   0.257   0.003   0.261
all.equal(taxa_names(ps1), taxa_names(ps2))
#> [1] TRUE
all.equal(otu_table(ps1), otu_table(ps2))
#> [1] TRUE
all.equal(tax_table(ps1), tax_table(ps2))
#> [1] TRUE
all.equal(phy_tree(ps1), phy_tree(ps2))
#> [1] TRUE
```

## Notes

My aim is for these functions to be drop-in replacements for phyloseq's
versions. See the unit tests in the `tests/` folder for the tests for
functional equivalence I've implemented so far. The `tax_glom()` functions
should produce identical results to phyloseq's version, though I have only
implemented tests using the default `NArm` and `bad_empty` arguments. The
`psmelt()` function in `phyloseq` drops columns of taxonomy data that are all
`NA` (such as after performing a `tax_glom()`), and returns a data frame with
extraneous row names. Speedyseq's `psmelt()` will not drop columns and does not
return row names. Both functions sort rows by the `Abundance` and `OTU`
columns, but the row order can differ in cases of ties for both variables.  If
you have any problems or find any discrepancies from `phyloseq`'s behavior,
please post an issue. Warning: Like phyloseq, speedyseq's `psmelt()` will
convert your taxonomy variables to factors if `getOption("stringsAsFactors")`
is `TRUE`.
