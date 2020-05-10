# speedyseq (development version)

* New `merge_taxa_vec()` function provides a vectorized version of
  `phyloseq::merge_taxa()` for `phyloseq` and `otu_table` objects.

* New `tip_glom()` function provides a speedy version of
  `phyloseq::tip_glom()` for indirect phylogenetic merging of taxa.

* New `tree_glom()` function performs direct phylogenetic merging of taxa. This
  function is much faster and arguably more intuitive than `tip_glom()`.

* Adds dependencies
  [castor](https://cran.r-project.org/web/packages/castor/index.html) and
  [purrr](https://purrr.tidyverse.org/)

## New general-purpose vectorized merging function

Phyloseq's `merge_taxa()` takes a phyloseq object or component object `x` and a
set of taxa `eqtaxa` and merges them into a single taxon. In place of the
`eqtaxa` argument, speedyseq's `merge_taxa_vec()` takes a vector `group` of
length `ntaxa(physeq)` that defines how all the taxa in `x` should be merged
into multiple new groups. Its syntax and behavior is patterned after that of
`base::rowsum()`, which it uses to do the merging in the OTU table. When aiming
to merge a large number of taxa into a smaller but still large number of
groups, it is much faster to do all the merging with one call to
`merge_taxa_vec()` than to loop through many calls to `merge_taxa()`.

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
tax_table(ps2)[c(108, 136, 45),]
```

## Speedy `tip_glom()` for indirect phylogenetic merging

Phyloseq provides `tip_glom()` to perform a form of indirect phylogenetic
merging using the phylogenetic tree in `phy_tree(physeq)`. This function uses
the tree to create a distance matrix, performs hierarchical clustering on the
distance matrix, and then defines new taxonomic groups by cutting the
dendrogram produced by the clustering at a user defined height. Phyloseq's
version can be slow and memory intensive when the number of taxa is large.

Speedyseq's new `tip_glom()` function provides a faster and less
memory-intensive alternative to `phyloseq::tip_glom() through the use of
vectorized merging (via `merge_taxa_vec()`) and faster and lower-memory
phylogenetic-distance computation (via `get_all_pairwise_distances()` from the
[castor](https://cran.r-project.org/web/packages/castor/index.html) package).

Speedyseq's `tax_glom()` also has the new `tax_adjust` argument, which is
passed on to `merge_taxa_vec()`. It is set to `1` by default for phyloseq
compatibility and should give identical results to phyloseq in this case.

For phyloseq compatibility, the default clustering function is left as
`cluster::agnes`. However, equivalent but faster results can be obtained by
using the `hclust` function from base R with the `method == "average"` option.

Speedyseq's `tip_glom()` currently only works on phyloseq objects and will give
an error if used on a phylo (tree) object.

## Direct phylogenetic merging with `tree_glom()`

It might be desirable in many cases to perform phylogenetic merging based
directly on the phylogenetic tree rather than (as in `tip_glom()`) a dendrogram
derived from it. Speedyseq's new `tree_glom()` function performs such direct
phylogenetic merging, which has several advantages.

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
