# speedyseq (development version)

* dplyr verbs for added for sample data and taxonomy tables: rename, rename_with, mutate, and filter (see #69)

* `merge_samples2()` now has a `fun_otu` argument for specifying alternative abundance-summarization functions

* Add  `ps_tibble()` S4 generic to provide `tibble::as_tibble()` functionality for phyloseq objects

* Fixed bug in print outputs. Row numbers are now kept as their removal was causing the issue. This is a temporary fix; see [#60](https://github.com/mikemc/speedyseq/issues/60).

* Fixed namespace bug that caused `psmelt(as="tibble")` to throw an error if
  tibble wasn't loaded

* `psmelt()` now uses `getOption(speedyseq.psmelt_class)` as the default value
  for the `as` argument. Users can set their preferred tabular class to 
  "speedyseq.psmelt_class" (among "data.table", "data.frame", or "tbl_df") in
  their ".Rprofile" file. The default option is set to "data.frame" for
  backwards compatibility.


# speedyseq 0.5.3

* `print(physeq)` now aligns the `refseq(physeq)` summary properly

# speedyseq 0.5.2

* New `as` argument in `psmelt()` allows specifying whether the result should
  be given as a "data.table", "data.frame", or "tbl_df" (tibble).

* In addition, `psmelt()` now ignores `options("stringsAsFactors")` and should
  never convert taxonomy to factors, in line with the phasing out of the
  `stringsAsFactors` behavior starting in R 4.0.0.

# speedyseq 0.5.1

* New function `orient_taxa()` to facilitate putting a phyloseq or otu-table
  object in a specific orientation (taxa as rows or as columns). This is useful
  when passing the otu table on to functions that require the abundance matrix
  to have a specific orientation and are unaware of the `taxa_are_rows(x)`
  property.

# speedyseq 0.5.0

* New "tibble-like" `show()` and `print()` methods for otu tables, sample-data
  tables, and taxonomy tables.

# speedyseq 0.4.1

* Fixed bug in `merge_samples2()` when the new sample names are the numerical
  sequence `1:n_groups`.

# speedyseq 0.4.0

* New `merge_samples2()` and helper `unique_or_na()` provides an alternative to
  `phyloseq::merge_samples()` that better handles categorical sample variables.
  The `funs` argument specifies which summary is used to merge each sample
  variable within groups. The default is `unique_or_na()`, which collapses the
  values to a single unique value if it exists and otherwise returns NA.

```r
data(enterotype)

# Merge samples with the same project and clinical status
ps <- enterotype 
sample_data(ps) <- sample_data(ps) %>%
  transform(Project.ClinicalStatus = Project:ClinicalStatus)
sample_data(ps) %>% head
ps0 <- merge_samples2(ps, "Project.ClinicalStatus", funs = list(Age = mean))
sample_data(ps0) %>% head
```

* Add whether taxa are rows to `show()` method for `phyloseq` objects

* Add `tibble::glimpse()` methods for `sample_data` and `phyloseq` objects

* Minor bug fixes to `merge_taxa_vec()`

# speedyseq 0.3.2

* Extend the constructor functions `otu_table()`, `sample_data()`, and
  `tax_table()` to work on [tibbles](https://r4ds.had.co.nz/tibbles.html).

# speedyseq 0.3.1

* Rename the function argument in `transform_sample_counts()` and
  `filter_taxa()` from `.f` to `fun` to match phyloseq.

# speedyseq 0.3.0

* New `transform_sample_counts()` and `filter_taxa()` provide wrappers around
  `phyloseq::transform_sample_counts()` and `phyloseq::filter_taxa()` that
  allow allow [purrr](https://purrr.tidyverse.org/)-style anonymous functions.

* New `filter_taxa2()` provides a version of `filter_taxa()` with `prune =
  TRUE`; that is, it always returns a pruned (filtered) phyloseq object. This
  version is convenient when filtering taxa in a pipe chain,

``` r
data(GlobalPatterns)
# Filter low prevalence taxa and then convert to proportions
gp.prop <- GlobalPatterns %>%
  filter_taxa2(~ sum(. > 0) > 5) %>%
  transform_sample_counts(~ . / sum(.))
```

* The [magrittr](https://magrittr.tidyverse.org/) pipe (`%>%`) is now exported
  so that it can be used without first loading magrittr or dplyr

# speedyseq 0.2.0

## Breaking changes

* The default ordering of new taxa output by `tax_glom()` is different from
  previous versions and from `phyloseq::tax_glom()` in phyloseq objects that do
  not have phylogenetic trees. See "Minor improvements and fixes" for more
  information.

## New features

### New general-purpose vectorized taxa-merging function

The new `merge_taxa_vec()` function provides a vectorized version of
`phyloseq::merge_taxa()` that can quickly merge arbitrary groups of taxa and
now forms the basis of all other merging functions. `phyloseq::merge_taxa()`
takes a phyloseq object or component object `x` and a set of taxa `eqtaxa` and
merges them into a single taxon. In place of the `eqtaxa` argument, speedyseq's
`merge_taxa_vec()` takes a vector `group` of length `ntaxa(physeq)` that
defines how all the taxa in `x` should be merged into multiple new groups. Its
syntax and behavior is patterned after that of `base::rowsum()`, which it uses
to do the merging in the OTU table. When aiming to merge a large number of taxa
into a smaller but still large number of groups, it is much faster to do all
the merging with one call to `merge_taxa_vec()` than to loop through many calls
to `merge_taxa()`.

A practical example is clustering amplicon sequence variants (ASVs) into OTUs
defined by a given similarity threshold. Suppose we have a phyloseq object `ps`
that has the ASV sequences stored in its `refseq` slot. We can cluster the ASV
sequences into 97% OTUs using the DECIPHER package with
```r
dna <- refseq(ps)
nproc <- 1 # Increase to use multiple processors
aln <- DECIPHER::AlignSeqs(dna, processors = nproc)
d <- DECIPHER::DistanceMatrix(aln, processors = nproc)
clusters <- DECIPHER::IdClusters(
  d, 
  method = "complete",
  cutoff = 0.03, # corresponds to 97% OTUs
  processors = nproc
)
```
Next we `merge_taxa_vec()` to get the merged phyloseq object,
```r
ps0 <- merge_taxa_vec(
  ps, 
  group = clusters$cluster,
  tax_adjust = 2
)
```
The names of the new taxa are set to the name of the most abundant taxon within
each group (the same behavior as the `tax_glom()` and `tip_glom()` functions).
Future versions will likely have a `names` argument to control the naming
behavior. 

The `tax_adjust` argument controls how NAs and within-group
disagreements in taxonomy are handled to determine the taxonomy of the new taxa
(see `help(merge_taxa_vec)` for details). An example of the difference between
`tax_adjust = 1` (`phyloseq::merge_taxa()` behavior) and `tax_adjust = 2` can
be seen in the following example from the new `tip_glom()` documentation,
```r
data(GlobalPatterns)

set.seed(20190421)
ps <- prune_taxa(sample(taxa_names(GlobalPatterns), 2e2), GlobalPatterns)

ps1 <- tip_glom(ps, 0.1, tax_adjust = 1)
ps2 <- tip_glom(ps, 0.1, tax_adjust = 2)
tax_table(ps1)[c(108, 136, 45),]
#> Taxonomy Table:     [3 taxa by 7 taxonomic ranks]:
#>        Kingdom    Phylum           Class                 Order               
#> 578831 "Bacteria" "Bacteroidetes"  "Sphingobacteria"     "Sphingobacteriales"
#> 2801   "Bacteria" "Planctomycetes" "Planctomycea"        "Pirellulales"      
#> 185581 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Oceanospirillales" 
#>        Family Genus            Species                             
#> 578831 NA     "Niabella"       NA                                  
#> 2801   NA     "Rhodopirellula" NA                                  
#> 185581 "OM60" NA               "marinegammaproteobacteriumHTCC2080"
tax_table(ps2)[c(108, 136, 45),]
#> Taxonomy Table:     [3 taxa by 7 taxonomic ranks]:
#>        Kingdom    Phylum           Class                 Order               
#> 578831 "Bacteria" "Bacteroidetes"  "Sphingobacteria"     "Sphingobacteriales"
#> 2801   "Bacteria" "Planctomycetes" "Planctomycea"        "Pirellulales"      
#> 185581 "Bacteria" "Proteobacteria" "Gammaproteobacteria" "Oceanospirillales" 
#>        Family Genus Species
#> 578831 NA     NA    NA     
#> 2801   NA     NA    NA     
#> 185581 "OM60" NA    NA     
```

### Faster and lower-memory implementation of `phyloseq::tip_glom()`

The new `tip_glom()` function provides a speedy version of
`phyloseq::tip_glom()`. This function performs a form of indirect phylogenetic
merging of taxa using the phylogenetic tree in `phy_tree(physeq)` by 1) using
the tree to create a distance matrix, 2) performing hierarchical clustering on
the distance matrix, and 3) defining new taxonomic groups by cutting the
dendrogram at the height specified by the `h` parameter. Speedyseq's
`tip_glom()` provides a faster and less memory-intensive alternative to
`phyloseq::tip_glom()` through the use of vectorized merging (via
`merge_taxa_vec()`) and faster and lower-memory phylogenetic-distance
computation (via `get_all_pairwise_distances()` from the
[castor](https://cran.r-project.org/web/packages/castor/index.html) package).

Speedyseq's `tip_glom()` also has the new `tax_adjust` argument, which is
passed on to `merge_taxa_vec()`. It is set to `1` by default for phyloseq
compatibility and should give identical results to phyloseq in this case.

For phyloseq compatibility, the default clustering function is left as
`cluster::agnes()`. However, equivalent but faster results can be obtained by
using the `hclust` function from base R with the `method == "average"` option.

### Direct phylogenetic merging with `tree_glom()`

The new `tree_glom()` function performs direct phylogenetic merging of taxa.
This function is much faster and arguably more intuitive than `tip_glom()`.
Advantages of direct merging over the indirect merging of `tip_glom()` are

1. A merged group of taxa correspond to a clade in the original tree being
   collapsed to a single taxon.
2. The `resolution` parameter that controls the degree of merging has
   units in terms of the tree's branch lengths, making it potentially more
   biologically meaningful than the `h` parameter in `tip_glom()`.
3. The distance-matrix computation and hierarchical clustering in `tip_glom()`
   can be skipped, making `tree_glom()` much faster and less memory intensive
   than `tip_glom()` when the number of taxa is large.

`tree_glom()` uses functions from the
[castor](https://cran.r-project.org/web/packages/castor/index.html) package to
determine which clades are to be merged using one of three criteria. The
default behavior is to merge a clade if the maximum distance from a node to
each of its tips is less than the distance `resolution`.

```r
data(GlobalPatterns)
ps1 <- subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
ntaxa(ps1)
ps2 <- tree_glom(ps1, 0.05)
ntaxa(ps2)

library(dplyr)
library(ggtree)
library(cowplot)

plot1 <- phy_tree(ps1) %>%
  ggtree +
  geom_tiplab() +
  geom_label(aes(x = branch, label = round(branch.length, 4)))
plot2 <- phy_tree(ps2) %>%
  ggtree +
  geom_tiplab() +
  geom_label(aes(x = branch, label = round(branch.length, 4)))
plot_grid(plot1, plot2)
```

## Minor improvements and fixes

* Fixed bug that applied to taxonomic merge functions when an object named
  `new_tax_mat` exists outside the function environment; described in 
  [Issue #31](https://github.com/mikemc/speedyseq/issues/31)

* Merging functions now maintain the original order of new taxa by default,
  except in phyloseq objects with phylogenetic trees (for which order is and
  has always been determined by how archetypes are ordered in
  `phy_tree(ps)$tip.label`). This behavior can lead to different taxa orders
  from past speedyseq versions and from `phyloseq::tax_glom()` function.
  However, it makes the resulting taxa order more predictable. New taxa can be
  be reordered according to `group` or taxonomy in `merge_taxa_vec()` and
  `tax_glom()` by setting `reorder = TRUE`.

* Merging/glom functions now work on relevant phyloseq components as well as
  phyloseq objects

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

`phyloseq::tax_glom()` can be applied to `taxonomyTable` objects as well as
`phyloseq` objects, but speedyseq's `tax_glom()` currently only works on
`phyloseq` objects and gives an error on `taxonomyTable` objects.
