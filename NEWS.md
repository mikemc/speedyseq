# speedyseq 0.1.2

* `tax_glom()` has a new implementation using base R functions instead of
  tibble and dplyr.

The new `tax_glom()` code takes advantage of the fact that we are working with
matrix data with a fixed row/col order and so do not need joining operations.
In particular, it uses `base::rowsum()` for simpler and faster merging taxa in
the OTU table (thanks to @digitalwright for the suggestion). Whether a
significant speed increase occurs depends on the dataset and call, since the
preliminary `prune_taxa()` step and the transpose step (to taxa-as-rows
orientation) are often the limiting steps and remain unchanged. Noticeable
speed ups on large phyloseq objects (~2x) can occur when `taxa_are_rows = TRUE`
(so that transposing is unnecessary).

# speedyseq 0.1.1

* `psmelt()` now uses data.table under the hood instead of the tidyverse
  packages tibble, tidyr, and dplyr. This refactor gives a roughly 2x speedup
  on the GlobalPatterns dataset.

# speedyseq 0.1.0

* New `psmelt()` uses functions from tidyr and dplyr in place of `merge()` and
  `reshape2::melt()` packages to achieve a dramatic speed up over
  `phyloseq::psmelt()` on large datasets

* New `tax_glom()` performs vectorized merging of taxonomic groups (using
  tibble and dplyr) to achieve a dramatic speed up over `phyloseq::tax_glom()`
  on large datasets

* Copies of the phyloseq plotting functions `plot_bar()`, `plot_heatmap()`, and
  `plot_tree()` are included that use the faster internal `psmelt()`

## Differences from phyloseq behavior

For most purposes, these functions should work as drop-in replacements for
phyloseq's versions, but there are a few differences to be aware of.

### `psmelt()`

The `psmelt()` function in `phyloseq` drops columns of taxonomy data that are
all `NA` (such as after performing a `tax_glom()`), and returns a data frame
with extraneous row names. Speedyseq's `psmelt()` will not drop columns and
does not return row names. Both functions sort rows by the `Abundance` and
`OTU` columns, but the row order can differ in cases of ties for both
variables. Warning: Like phyloseq's version, speedyseq's `psmelt()` will
convert your taxonomy variables to factors if `getOption("stringsAsFactors")`
is `TRUE`.

### `tax_glom()`

Phyloseq's `tax_glom()` can be applied to `taxonomyTable` objects as well as
`phyloseq` objects, but speedyseq's `tax_glom()` currently only works on
`phyloseq` objects and gives an error on `taxonomyTable` objects.
