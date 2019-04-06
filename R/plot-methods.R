# Note, these are all copied straight from phyloseq, with added importFrom's
# and calling of phyloseq internal functions as needed to adjust for the
# missing functions from phyloseq's plot-methods.R
################################################################################
#' A flexible, informative barplot phyloseq data
#'
#' There are many useful examples of phyloseq barplot graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_bar-examples}{phyloseq online tutorials}.
#' This function wraps \code{ggplot2} plotting, and returns a \code{ggplot2}
#'  graphic object
#' that can be saved or further modified with additional layers, options, etc.
#' The main purpose of this function is to quickly and easily create informative
#' summary graphics of the differences in taxa abundance between samples in
#' an experiment. 
#'
#' @usage plot_bar(physeq, x="Sample", y="Abundance", fill=NULL,
#'  title=NULL, facet_grid=NULL)
#'
#' @param physeq (Required). An \code{\link{otu_table-class}} or 
#'  \code{\link{phyloseq-class}}.
#'
#' @param x (Optional). Optional, but recommended, especially if your data
#'  is comprised of many samples. A character string.
#'  The variable in the melted-data that should be mapped to the x-axis.
#'  See \code{\link{psmelt}}, \code{\link{melt}},
#'  and \code{\link{ggplot}} for more details.
#' 
#' @param y (Optional). A character string.
#'  The variable in the melted-data that should be mapped to the y-axis.
#'  Typically this will be \code{"Abundance"}, in order to
#'  quantitatively display the abundance values for each OTU/group. 
#'  However, alternative variables could be used instead,
#'  producing a very different, though possibly still informative, plot.
#'  See \code{\link{psmelt}}, \code{\link{melt}},
#'  and \code{\link{ggplot}} for more details.
#'
#' @param fill (Optional). A character string. Indicates which sample variable
#'  should be used to map to the fill color of the bars. 
#'  The default is \code{NULL}, resulting in a gray fill for all bar segments.
#' 
#' @param facet_grid (Optional). A formula object.
#'  It should describe the faceting you want in exactly the same way as for 
#'  \code{\link[ggplot2]{facet_grid}}, 
#'  and is ulitmately provided to \code{\link{ggplot}}2 graphics.
#'  The default is: \code{NULL}, resulting in no faceting.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @return A \code{\link[ggplot2]{ggplot}}2 graphic object -- rendered in the graphical device
#'  as the default \code{\link[base]{print}}/\code{\link[methods]{show}} method.
#'
#' @seealso 
#'  \href{http://joey711.github.io/phyloseq/plot_bar-examples}{phyloseq online tutorials}.
#'
#'  \code{\link{psmelt}}
#'
#'  \code{\link{ggplot}}
#' 
#'  \code{\link{qplot}}
#'
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes_string
#' @importFrom ggplot2 geom_bar
#' @importFrom ggplot2 facet_grid
#' @importFrom ggplot2 element_text
#' @importFrom ggplot2 theme
#' @importFrom ggplot2 ggtitle
#' 
#' @export
#'
#' @examples
#' data("GlobalPatterns")
#' gp.ch = subset_taxa(GlobalPatterns, Phylum == "Chlamydiae")
#' plot_bar(gp.ch)
#' plot_bar(gp.ch, fill="Genus")
#' plot_bar(gp.ch, x="SampleType", fill="Genus")
#' plot_bar(gp.ch, "SampleType", fill="Genus", facet_grid=~Family)
#' # See additional examples in the plot_bar online tutorial. Link above.
plot_bar = function(physeq, x="Sample", y="Abundance", fill=NULL,
	title=NULL, facet_grid=NULL){
		
	# Start by melting the data in the "standard" way using psmelt.
	mdf = psmelt(physeq)
	
	# Build the plot data structure
	p = ggplot(mdf, aes_string(x=x, y=y, fill=fill))

	# Add the bar geometric object. Creates a basic graphic. Basis for the rest.
	# Test weather additional
	p = p + geom_bar(stat="identity", position="stack", color="black")

	# By default, rotate the x-axis labels (they might be long)
	p = p + theme(axis.text.x=element_text(angle=-90, hjust=0))

	# Add faceting, if given
	if( !is.null(facet_grid) ){	
		p <- p + facet_grid(facet_grid)
	}
	
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p <- p + ggtitle(title)
	}
	
	return(p)
}
################################################################################
#' Plot a phylogenetic tree with optional annotations
#'
#' There are many useful examples of phyloseq tree graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_tree-examples}{phyloseq online tutorials}.
#' This function is intended to facilitate easy graphical investigation of 
#' the phylogenetic tree, as well as sample data. Note that for phylogenetic
#' sequencing of samples with large richness, some of the options in this 
#' function will be prohibitively slow to render, or too dense to be
#' interpretable. A rough ``rule of thumb'' is to use subsets of data 
#' with not many more than 200 OTUs per plot, sometimes less depending on the
#' complexity of the additional annotations being mapped to the tree. It is 
#' usually possible to create an unreadable, uninterpretable tree with modern
#' datasets. However, the goal should be toward parameter settings and data
#' subsets that convey (honestly, accurately) some biologically relevant
#' feature of the data. One of the goals of the \code{\link{phyloseq-package}}
#' is to make the determination of these features/settings as easy as possible.
#'
#' This function received an early development contribution from the work of 
#' Gregory Jordan via \href{https://github.com/gjuggler/ggphylo}{the ggphylo package}.
#' \code{plot_tree} has since been re-written.
#' For details see \code{\link{tree_layout}}.
#'
#' @param physeq (Required). The data about which you want to 
#'  plot and annotate a phylogenetic tree, in the form of a
#'  single instance of the \code{\link{phyloseq-class}}, containing at 
#'  minimum a phylogenetic tree component (try \code{\link{phy_tree}}).
#'  One of the major advantages of this function over basic tree-plotting utilities
#'  in the \code{\link{ape}}-package is the ability to easily annotate the tree
#'  with sample variables and taxonomic information. For these uses, 
#'  the \code{physeq} argument should also have a \code{\link{sample_data}}
#'  and/or \code{\link{tax_table}} component(s).
#' 
#' @param method (Optional). Character string. Default \code{"sampledodge"}. 
#'  The name of the annotation method to use. 
#'  This will be expanded in future versions.
#'  Currently only \code{"sampledodge"} and \code{"treeonly"} are supported.
#'  The \code{"sampledodge"} option results in points
#'  drawn next to leaves if individuals from that taxa were observed,
#'  and a separate point is drawn for each sample.
#' 
#' @param nodelabf (Optional). A function. Default \code{NULL}.
#'  If \code{NULL}, the default, a function will be selected for you based upon
#'  whether or not there are node labels in \code{phy_tree(physeq)}.
#'  For convenience, the phyloseq package includes two generator functions
#'  for adding arbitrary node labels (can be any character string),
#'  \code{\link{nodeplotdefault}};
#'  as well as for adding bootstrap values in a certain range,
#'  \code{\link{nodeplotboot}}.
#'  To not have any node labels in the graphic, set this argument to
#'  \code{\link{nodeplotblank}}.
#'
#' @param color (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point color.
#'  Supported options here also include the reserved special variables
#'  of \code{\link{psmelt}}.
#' 
#' @param shape (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point shape.
#'  Supported options here also include the reserved special variables
#'  of \code{\link{psmelt}}.
#'
#' @param size (Optional). Character string. Default \code{NULL}.
#'  The name of the variable in \code{physeq} to map to point size.
#'  A special argument \code{"abundance"} is reserved here and scales
#'  point size using abundance in each sample on a log scale.
#'  Supported options here also include the reserved special variables
#'  of \code{\link{psmelt}}.
#'
#' @param min.abundance (Optional). Numeric. 
#'  The minimum number of individuals required to label a point
#'  with the precise number.
#'  Default is \code{Inf},
#'  meaning that no points will have their abundance labeled.
#'  If a vector, only the first element is used. 
#'
#' @param label.tips (Optional). Character string. Default is \code{NULL},
#'  indicating that no tip labels will be printed.
#'  If \code{"taxa_names"}, then the name of the taxa will be added 
#'  to the tree; either next to the leaves, or next to
#'  the set of points that label the leaves. Alternatively,
#'  if this is one of the rank names (from \code{rank_names(physeq)}),
#'  then the identity (if any) for that particular taxonomic rank
#'  is printed instead.
#'
#' @param text.size (Optional). Numeric. Should be positive. The 
#'  size parameter used to control the text size of taxa labels.
#'  Default is \code{NULL}. If left \code{NULL}, this function
#'  will automatically calculate a (hopefully) optimal text size
#'  given the vertical constraints posed by the tree itself. 
#'  This argument is included in case the 
#'  automatically-calculated size is wrong, and you want to change it.
#'  Note that this parameter is only meaningful if \code{label.tips}
#'  is not \code{NULL}.
#'
#' @param sizebase (Optional). Numeric. Should be positive.
#'  The base of the logarithm used
#'  to scale point sizes to graphically represent abundance of
#'  species in a given sample. Default is 5.
#' 
#' @param base.spacing (Optional). Numeric. Default is \code{0.02}.
#'  Should be positive.
#'  This defines the base-spacing between points at each tip/leaf in the
#'  the tree. The larger this value, the larger the spacing between points.
#'  This is useful if you have problems with overlapping large points
#'  and/or text indicating abundance, for example. Similarly, if you 
#'  don't have this problem and want tighter point-spacing, you can 
#'  shrink this value.
#'
#' @param ladderize (Optional). Boolean or character string (either
#'  \code{FALSE}, \code{TRUE}, or \code{"left"}).
#'  Default is \code{FALSE}.
#'  This parameter specifies whether or not to \code{\link[ape]{ladderize}} the tree 
#'  (i.e., reorder nodes according to the depth of their enclosed
#'  subtrees) prior to plotting.
#'  This tends to make trees more aesthetically pleasing and legible in
#'  a graphical display.
#'  When \code{TRUE} or \code{"right"}, ``right'' ladderization is used.
#'  When set to \code{FALSE}, no ladderization is applied.
#'  When set to \code{"left"}, the reverse direction
#'  (``left'' ladderization) is applied.
#'  This argument is passed on to \code{\link{tree_layout}}.
#'
#' @param plot.margin (Optional). Numeric. Default is \code{0.2}.
#'  Should be positive.
#'  This defines how much right-hand padding to add to the tree plot,
#'  which can be required to not truncate tip labels. The margin value
#'  is specified as a fraction of the overall tree width which is added
#'  to the right side of the plot area. So a value of \code{0.2} adds
#'  twenty percent extra space to the right-hand side of the plot.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'  
#' @param treetheme (Optional).
#'  A custom \code{\link{ggplot}}2 \code{\link[ggplot2]{theme}} layer
#'  to use for the tree. Supplants any default theme layers 
#'  used within the function.
#'  A value of \code{NULL} uses a default, minimal-annotations theme. 
#'  If anything other than a them or \code{NULL}, the current global ggplot2
#'  theme will result.
#'  
#' @param justify (Optional). A character string indicating the
#'  type of justification to use on dodged points and tip labels. 
#'  A value of \code{"jagged"}, the default, results in 
#'  these tip-mapped elements being spaced as close to the tips as possible
#'  without gaps. 
#'  Currently, any other value for \code{justify} results in
#'  a left-justified arrangement of both labels and points.
#'
#' @return A \code{\link{ggplot}}2 plot.
#' 
#' @seealso
#'  \code{\link{plot.phylo}}
#'
#' There are many useful examples of phyloseq tree graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_tree-examples}{phyloseq online tutorials}.
#'
#' @importFrom scales log_trans
#' 
#' @importFrom data.table data.table
#' @importFrom data.table setkey
#' @importFrom data.table setkeyv
#' 
#' @importFrom ggplot2 geom_segment
#' @importFrom ggplot2 geom_point
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 scale_size_continuous
#' @importFrom ggplot2 element_blank
#' 
#' @export
#' @examples
#' # # Using plot_tree() with the esophagus dataset.
#' # # Please note that many more interesting examples are shown
#' # # in the online tutorials"
#' # # http://joey711.github.io/phyloseq/plot_tree-examples
#' data(esophagus)
#' # plot_tree(esophagus)
#' # plot_tree(esophagus, color="Sample")
#' # plot_tree(esophagus, size="Abundance")
#' # plot_tree(esophagus, size="Abundance", color="samples")
#' plot_tree(esophagus, size="Abundance", color="Sample", base.spacing=0.03)
#' plot_tree(esophagus, size="abundance", color="samples", base.spacing=0.03)
plot_tree = function(physeq, method="sampledodge", nodelabf=NULL,
                       color=NULL, shape=NULL, size=NULL,
                       min.abundance=Inf, label.tips=NULL, text.size=NULL,
                       sizebase=5, base.spacing = 0.02,
                       ladderize=FALSE, plot.margin=0.2, title=NULL,
                       treetheme=NULL, justify="jagged"){
  ########################################
  # Support mis-capitalization of reserved variable names in color, shape, size
  # This helps, for instance, with backward-compatibility where "abundance"
  # was the reserved variable name for mapping OTU abundance entries
  fix_reserved_vars = function(aesvar){
    aesvar <- gsub("^abundance[s]{0,}$", "Abundance", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^OTU[s]{0,}$", "OTU", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^taxa_name[s]{0,}$", "OTU", aesvar, ignore.case=TRUE)
    aesvar <- gsub("^sample[s]{0,}$", "Sample", aesvar, ignore.case=TRUE)
    return(aesvar)
  }
  if(!is.null(label.tips)){label.tips <- fix_reserved_vars(label.tips)}
  if(!is.null(color)){color <- fix_reserved_vars(color)}
  if(!is.null(shape)){shape <- fix_reserved_vars(shape)}
  if(!is.null(size) ){size  <- fix_reserved_vars(size)} 
  ########################################
  if( is.null(phy_tree(physeq, FALSE)) ){
    stop("There is no phylogenetic tree in the object you have provided.\n",
         "Try phy_tree(physeq) to see for yourself.")
  }
  if(!inherits(physeq, "phyloseq")){
    # If only a phylogenetic tree, then only tree available to overlay.
    method <- "treeonly"
  }
  # Create the tree data.table
  treeSegs <- phyloseq::tree_layout(phy_tree(physeq), ladderize=ladderize)
  edgeMap = aes(x=xleft, xend=xright, y=y, yend=y)
  vertMap = aes(x=x, xend=x, y=vmin, yend=vmax)
  # Initialize phylogenetic tree.
  # Naked, lines-only, unannotated tree as first layers. Edge (horiz) first, then vertical.
  p = ggplot(data=treeSegs$edgeDT) + geom_segment(edgeMap) + 
    geom_segment(vertMap, data=treeSegs$vertDT)
  # If no text.size given, calculate it from number of tips ("species", aka taxa)
  # This is very fast. No need to worry about whether text is printed or not.
  if(is.null(text.size)){
    text.size <- phyloseq:::manytextsize(ntaxa(physeq))
  }
  # Add the species labels to the right.
  if(!is.null(label.tips) & method!="sampledodge"){
    # If method is sampledodge, then labels are added to the right of points, later.
    # Add labels layer to plotting object.
    labelDT = treeSegs$edgeDT[!is.na(OTU), ]
    if(!is.null(tax_table(object=physeq, errorIfNULL=FALSE))){
      # If there is a taxonomy available, merge it with the label data.table
      taxDT = data.table(tax_table(physeq), OTU=taxa_names(physeq), key="OTU")
      # Merge with taxonomy.
      labelDT = merge(x=labelDT, y=taxDT, by="OTU")
    }
    if(justify=="jagged"){
      # Tip label aesthetic mapping.
      # Aesthetics can be NULL, and that aesthetic gets ignored.
      labelMap <- aes_string(x="xright", y="y", label=label.tips, color=color)
    } else {
      # The left-justified version of tip-labels.
      labelMap <- aes_string(x="max(xright, na.rm=TRUE)", y="y", label=label.tips, color=color)
    }
    p <- p + geom_text(labelMap, data=labelDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
  }
  # Node label section.
  # 
  # If no nodelabf ("node label function") given, ask internal function to pick one.
  # Is NULL by default, meaning will dispatch to `howtolabnodes` to select function.
  # For no node labels, the "dummy" function `nodeplotblank` will return tree plot 
  # object, p, as-is, unmodified.
  if(is.null(nodelabf)){
    nodelabf = phyloseq:::howtolabnodes(physeq)
  }
  #### set node `y` as the mean of the vertical segment
  # Use the provided/inferred node label function to add the node labels layer(s)
  # Non-root nodes first
  p = nodelabf(p, treeSegs$edgeDT[!is.na(label), ])
  # Add root label (if present)
  p = nodelabf(p, treeSegs$vertDT[!is.na(label), ])
  # Theme specification
  if(is.null(treetheme)){
    # If NULL, then use the default tree theme.
    treetheme <- theme(axis.ticks = element_blank(),
                       axis.title.x=element_blank(), axis.text.x=element_blank(),
                       axis.title.y=element_blank(), axis.text.y=element_blank(),
                       panel.background = element_blank(),
                       panel.grid.minor = element_blank(),      
                       panel.grid.major = element_blank())   
  }
  if(inherits(treetheme, "theme")){
    # If a theme, add theme layer to plot. 
    # For all other cases, skip this, which will cause default theme to be used
    p <- p + treetheme
  }
  # Optionally add a title to the plot
  if(!is.null(title)){
    p <- p + ggtitle(title)
  }  
  if(method!="sampledodge"){
    # If anything but a sampledodge tree, return now without further decorations.
    return(p)
  }
  ########################################
  # Sample Dodge Section
  # Special words, c("Sample", "Abundance", "OTU")
  # See psmelt()
  ########################################
  # Initialize the species/taxa/OTU data.table
  dodgeDT = treeSegs$edgeDT[!is.na(OTU), ]
  # Merge with psmelt() result, to make all co-variables available
  dodgeDT = merge(x=dodgeDT, y=data.table(psmelt(physeq), key="OTU"), by="OTU")
  if(justify=="jagged"){
    # Remove 0 Abundance value entries now, not later, for jagged.
    dodgeDT <- dodgeDT[Abundance > 0, ]    
  }
  # Set key. Changes `dodgeDT` in place. OTU is first key, always.
  if( !is.null(color) | !is.null(shape) | !is.null(size) ){
    # If color, shape, or size is chosen, setkey by those as well
    setkeyv(dodgeDT, cols=c("OTU", color, shape, size))
  } else {
    # Else, set key by OTU and sample name. 
    setkey(dodgeDT, OTU, Sample)
  }
  # Add sample-dodge horizontal adjustment index. In-place data.table assignment
  dodgeDT[, h.adj.index := 1:length(xright), by=OTU]
  # `base.spacing` is a user-input parameter.
  # The sampledodge step size is based on this and the max `x` value
  if(justify=="jagged"){
    dodgeDT[, xdodge:=(xright + h.adj.index * base.spacing * max(xright, na.rm=TRUE))]
  } else {
    # Left-justified version, `xdodge` always starts at the max of all `xright` values.
    dodgeDT[, xdodge := max(xright, na.rm=TRUE) + h.adj.index * base.spacing * max(xright, na.rm=TRUE)]
    # zeroes removed only after all sample points have been mapped.
    dodgeDT <- dodgeDT[Abundance > 0, ]
  }
  # The general tip-point map. Objects can be NULL, and that aesthetic gets ignored.
  dodgeMap <- aes_string(x="xdodge", y="y", color=color, fill=color,
                        shape=shape, size=size)
  p <- p + geom_point(dodgeMap, data=dodgeDT, na.rm=TRUE)
  # Adjust point size transform
  if( !is.null(size) ){
    p <- p + scale_size_continuous(trans=log_trans(sizebase))
  }  
  # Optionally-add abundance value label to each point.
  # User controls this by the `min.abundance` parameter.
  # A value of `Inf` implies no labels.
  if( any(dodgeDT$Abundance >= min.abundance[1]) ){
    pointlabdf = dodgeDT[Abundance>=min.abundance[1], ]
    p <- p + geom_text(mapping=aes(xdodge, y, label=Abundance),
                       data=pointlabdf, size=text.size, na.rm=TRUE)
  }
  # If indicated, add the species labels to the right of dodged points.
  if(!is.null(label.tips)){
    # `tiplabDT` has only one row per tip, the farthest horizontal
    # adjusted position (one for each taxa)
    tiplabDT = dodgeDT
    tiplabDT[, xfartiplab:=max(xdodge), by=OTU]
    tiplabDT <- tiplabDT[h.adj.index==1, .SD, by=OTU]
    if(!is.null(color)){
      if(color %in% sample_variables(physeq, errorIfNULL=FALSE)){
        color <- NULL
      }
    }
    labelMap <- NULL
    if(justify=="jagged"){
      labelMap <- aes_string(x="xfartiplab", y="y", label=label.tips, color=color)
    } else {
      labelMap <- aes_string(x="max(xfartiplab, na.rm=TRUE)", y="y", label=label.tips, color=color)
    }
    # Add labels layer to plotting object.
    p <- p + geom_text(labelMap, tiplabDT, size=I(text.size), hjust=-0.1, na.rm=TRUE)
  } 
  # Plot margins. 
  # Adjust the tree graphic plot margins.
  # Helps to manually ensure that graphic elements aren't clipped,
  # especially when there are long tip labels.
  min.x <- -0.01 # + dodgeDT[, min(c(xleft))]
  max.x <- dodgeDT[, max(xright, na.rm=TRUE)]
  if("xdodge" %in% names(dodgeDT)){
    max.x <- dodgeDT[, max(xright, xdodge, na.rm=TRUE)]
  }
  if(plot.margin > 0){
    max.x <- max.x * (1.0 + plot.margin)
  } 
  p <- p + scale_x_continuous(limits=c(min.x, max.x))  
  return(p)
}
################################################################################
#' Create an ecologically-organized heatmap using ggplot2 graphics
#'
#' There are many useful examples of phyloseq heatmap graphics in the
#' \href{http://joey711.github.io/phyloseq/plot_heatmap-examples}{phyloseq online tutorials}.
#' In a 2010 article in BMC Genomics, Rajaram and Oono show describe an 
#' approach to creating a heatmap using ordination methods to organize the 
#' rows and columns instead of (hierarchical) cluster analysis. In many cases
#' the ordination-based ordering does a much better job than h-clustering. 
#' An immediately useful example of their approach is provided in the NeatMap
#' package for R. The NeatMap package can be used directly on the abundance 
#' table (\code{\link{otu_table-class}}) of phylogenetic-sequencing data, but 
#' the NMDS or PCA ordination options that it supports are not based on ecological
#' distances. To fill this void, phyloseq provides the \code{plot_heatmap()}
#' function as an ecology-oriented variant of the NeatMap approach to organizing
#' a heatmap and build it using ggplot2 graphics tools.
#' The \code{distance} and \code{method} arguments are the same as for the
#' \code{\link{plot_ordination}} function, and support large number of
#' distances and ordination methods, respectively, with a strong leaning toward
#' ecology.
#' This function also provides the options to re-label the OTU and sample 
#' axis-ticks with a taxonomic name and/or sample variable, respectively, 
#' in the hope that this might hasten your interpretation of the patterns
#' (See the \code{sample.label} and \code{taxa.label} documentation, below). 
#' Note that this function makes no attempt to overlay hierarchical 
#' clustering trees on the axes, as hierarchical clustering is not used to 
#' organize the plot. Also note that each re-ordered axis repeats at the edge,
#' and so apparent clusters at the far right/left or top/bottom of the 
#' heat-map may actually be the same. For now, the placement of this edge
#' can be considered arbitrary, so beware of this artifact of this graphical
#' representation. If you benefit from this phyloseq-specific implementation
#' of the NeatMap approach, please cite both our packages/articles.
#'
#' This approach borrows heavily from the \code{heatmap1} function in the
#' \code{NeatMap} package. Highly recommended, and we are grateful for their
#' package and ideas, which we have adapted for our specific purposes here,
#' but did not use an explicit dependency. At the time of the first version
#' of this implementation, the NeatMap package depends on the rgl-package,
#' which is not needed in phyloseq, at present. Although likely a transient
#' issue, the rgl-package has some known installation issues that have further
#' influenced to avoid making NeatMap a formal dependency (Although we love
#' both NeatMap and rgl!).
#'
#' @param physeq (Required). The data, in the form of an instance of the
#'  \code{\link{phyloseq-class}}. This should be what you get as a result
#'  from one of the
#'  \code{\link{import}} functions, or any of the processing downstream.
#'  No data components beyond the \code{\link{otu_table}} are strictly 
#'  necessary, though they may be useful if you want to re-label the 
#'  axis ticks according to some observable or taxonomic rank, for instance,
#'  or if you want to use a \code{\link{UniFrac}}-based distance
#'  (in which case your \code{physeq} data would need to have a tree included).
#'
#' @param method (Optional).
#'  The ordination method to use for organizing the 
#'  heatmap. A great deal of the usefulness of a heatmap graphic depends upon 
#'  the way in which the rows and columns are ordered. 
#'
#' @param distance (Optional). A character string. 
#'  The ecological distance method to use in the ordination.
#'  See \code{\link{distance}}.
#'
#' @param sample.label (Optional). A character string.
#'  The sample variable by which you want to re-label the sample (horizontal) axis.
#'
#' @param taxa.label (Optional). A character string.
#'  The name of the taxonomic rank by which you want to
#'  re-label the taxa/species/OTU (vertical) axis.
#'  You can see available options in your data using
#'  \code{\link{rank_names}(physeq)}.
#'
#' @param low (Optional). A character string. An R color.
#'  See \code{?\link{colors}} for options support in R (there are lots).
#'  The color that represents the lowest non-zero value
#'  in the heatmap. Default is a dark blue color, \code{"#000033"}.
#' 
#' @param high (Optional). A character string. An R color.
#'  See \code{\link{colors}} for options support in R (there are lots).
#'  The color that will represent the highest 
#'  value in the heatmap. The default is \code{"#66CCFF"}.
#'  Zero-values are treated as \code{NA}, and set to \code{"black"}, to represent
#'  a background color.
#'
#' @param na.value (Optional). A character string. An R color.
#'  See \code{\link{colors}} for options support in R (there are lots).
#'  The color to represent what is essentially the background of the plot,
#'  the non-observations that occur as \code{NA} or
#'  \code{0} values in the abundance table. The default is \code{"black"}, which 
#'  works well on computer-screen graphics devices, but may be a poor choice for
#'  printers, in which case you might want this value to be \code{"white"}, and
#'  reverse the values of \code{high} and \code{low}, above.
#'
#' @param trans (Optional). \code{"trans"}-class transformer-definition object.
#'  A numerical transformer to use in 
#'  the continuous color scale. See \code{\link[scales]{trans_new}} for details.
#'  The default is \code{\link{log_trans}(4)}.
#'
#' @param max.label (Optional). Integer. Default is 250.
#'  The maximum number of labeles to fit on a given axis (either x or y). 
#'  If number of taxa or samples exceeds this value, 
#'  the corresponding axis will be stripped of any labels. 
#'
#'  This supercedes any arguments provided to
#'  \code{sample.label} or \code{taxa.label}.
#'  Make sure to increase this value if, for example,
#'  you want a special label
#'  for an axis that has 300 indices.
#'
#' @param title (Optional). Default \code{NULL}. Character string.
#'  The main title for the graphic.
#'
#' @param sample.order (Optional). Default \code{NULL}. 
#'  Either a single character string matching 
#'  one of the \code{\link{sample_variables}} in your data,
#'  or a character vector of \code{\link{sample_names}}
#'  in the precise order that you want them displayed in the heatmap.
#'  This overrides any ordination ordering that might be done
#'  with the \code{method}/\code{distance} arguments.
#'  
#' @param taxa.order (Optional). Default \code{NULL}. 
#'  Either a single character string matching 
#'  one of the \code{\link{rank_names}} in your data,
#'  or a character vector of \code{\link{taxa_names}}
#'  in the precise order that you want them displayed in the heatmap.
#'  This overrides any ordination ordering that might be done
#'  with the \code{method}/\code{distance} arguments.
#' 
#' @param first.sample (Optional). Default \code{NULL}.
#'  A character string matching one of the \code{\link{sample_names}}
#'  of your input data (\code{physeq}). 
#'  It will become the left-most sample in the plot.
#'  For the ordination-based ordering (recommended),
#'  the left and right edges of the axes are adjaacent in a continuous ordering. 
#'  Therefore, the choice of starting sample is meaningless and arbitrary,
#'  but it is aesthetically poor to have the left and right edge split 
#'  a natural cluster in the data.
#'  This argument allows you to specify the left edge
#'  and thereby avoid cluster-splitting, emphasize a gradient, etc.
#'  
#' @param first.taxa (Optional). Default \code{NULL}.
#'  A character string matching one of the \code{\link{taxa_names}}
#'  of your input data (\code{physeq}). 
#'  This is equivalent to \code{first.sample} (above),
#'  but for the taxa/OTU indices, usually the vertical axis.
#'
#' @param ... (Optional). Additional parameters passed to \code{\link{ordinate}}.
#' 
#' @return
#'  A heatmap plot, in the form of a \code{\link{ggplot}2} plot object,
#'  which can be further saved and modified.
#'
#' @references
#'  Because this function relies so heavily in principle, and in code, on some of the
#'  functionality in NeatMap, please site their article if you use this function
#'  in your work.
#' 
#'  Rajaram, S., & Oono, Y. (2010).
#'  NeatMap--non-clustering heat map alternatives in R. BMC Bioinformatics, 11, 45.
#'
#' Please see further examples in the 
#' \href{http://joey711.github.io/phyloseq/plot_heatmap-examples}{phyloseq online tutorials}.
#' 
#' @importFrom vegan scores
#' 
#' @importFrom scales log_trans
#' 
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 scale_y_discrete
#' @importFrom ggplot2 scale_x_discrete
#' @importFrom ggplot2 scale_fill_gradient
#' @importFrom ggplot2 geom_raster
#' @importFrom ggplot2 theme
#' 
#' @export
#' @examples
#' data("GlobalPatterns")
#' gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
#' # FYI, the base-R function uses a non-ecological ordering scheme,
#' # but does add potentially useful hclust dendrogram to the sides...
#' gpac <- subset_taxa(GlobalPatterns, Phylum=="Crenarchaeota")
#' # Remove the nearly-empty samples (e.g. 10 reads or less)
#' gpac = prune_samples(sample_sums(gpac) > 50, gpac)
#' # Arbitrary order if method set to NULL
#' plot_heatmap(gpac, method=NULL, sample.label="SampleType", taxa.label="Family")
#' # Use ordination
#' plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family")
#' # Use ordination for OTUs, but not sample-order
#' plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family", sample.order="SampleType")
#' # Specifying both orders omits any attempt to use ordination. The following should be the same.
#' p0 = plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family", taxa.order="Phylum", sample.order="SampleType")
#' p1 = plot_heatmap(gpac, method=NULL, sample.label="SampleType", taxa.label="Family", taxa.order="Phylum", sample.order="SampleType")
#' #expect_equivalent(p0, p1)
#' # Example: Order matters. Random ordering of OTU indices is difficult to interpret, even with structured sample order
#' rando = sample(taxa_names(gpac), size=ntaxa(gpac), replace=FALSE)
#' plot_heatmap(gpac, method=NULL, sample.label="SampleType", taxa.label="Family", taxa.order=rando, sample.order="SampleType")
#' # # Select the edges of each axis. 
#' # First, arbitrary edge, ordering
#' plot_heatmap(gpac, method=NULL)
#' # Second, biological-ordering (instead of default ordination-ordering), but arbitrary edge
#' plot_heatmap(gpac, taxa.order="Family", sample.order="SampleType")
#' # Third, biological ordering, selected edges
#' plot_heatmap(gpac, taxa.order="Family", sample.order="SampleType", first.taxa="546313", first.sample="NP2")
#' # Fourth, add meaningful labels
#' plot_heatmap(gpac, sample.label="SampleType", taxa.label="Family", taxa.order="Family", sample.order="SampleType", first.taxa="546313", first.sample="NP2")
plot_heatmap <- function(physeq, method="NMDS", distance="bray", 
	sample.label=NULL, taxa.label=NULL, 
	low="#000033", high="#66CCFF", na.value="black", trans=log_trans(4), 
	max.label=250, title=NULL, sample.order=NULL, taxa.order=NULL,
  first.sample=NULL, first.taxa=NULL, ...){

  # User-override ordering
  if( !is.null(taxa.order) & length(taxa.order)==1 ){
    # Assume `taxa.order` is a tax_table variable. Use it for ordering.
    rankcol = which(rank_names(physeq) %in% taxa.order)
    taxmat = as(tax_table(physeq)[, 1:rankcol], "matrix")
    taxa.order = apply(taxmat, 1, paste, sep="", collapse="")
    names(taxa.order) <- taxa_names(physeq)
    taxa.order = names(sort(taxa.order, na.last=TRUE))
  }
  if( !is.null(sample.order) & length(sample.order)==1 ){
    # Assume `sample.order` is a sample variable. Use it for ordering.
    sample.order = as.character(get_variable(physeq, sample.order))
    names(sample.order) <- sample_names(physeq)
    sample.order = names(sort(sample.order, na.last=TRUE))
  }
  
  if( !is.null(method) & (is.null(taxa.order) | is.null(sample.order)) ){
    # Only attempt NeatMap if method is non-NULL & at least one of
    # taxa.order and sample.order is not-yet defined.
    # If both axes orders pre-defined by user, no need to perform ordination...
    
	  # Copy the approach from NeatMap by doing ordination on samples, but use 
	  # phyloseq-wrapped distance/ordination procedures.
	  # Reorder by the angle in radial coordinates on the 2-axis plane.
    
    # In case of NMDS iterations, capture the output so it isn't dumped on standard-out
		junk = capture.output( ps.ord <- ordinate(physeq, method, distance, ...), file=NULL)
    if( is.null(sample.order) ){
      siteDF = NULL
      # Only define new ord-based sample order if user did not define one already
      trash1 = try({siteDF <- scores(ps.ord, choices = c(1, 2), display="sites", physeq=physeq)},
                   silent = TRUE)
      if(inherits(trash1, "try-error")){
        # warn that the attempt to get ordination coordinates for ordering failed.
        warning("Attempt to access ordination coordinates for sample ordering failed.\n",
                "Using default sample ordering.")
      }
      if(!is.null(siteDF)){
        # If the score accession seemed to work, go ahead and replace sample.order
        sample.order <- sample_names(physeq)[order(phyloseq:::RadialTheta(siteDF))]
      }
    }

		if( is.null(taxa.order) ){
		  # re-order species/taxa/OTUs, if possible,
		  # and only if user did not define an order already
			specDF = NULL
			trash2 = try({specDF <- scores(ps.ord, choices=c(1, 2), display="species", physeq=physeq)},
			             silent = TRUE)
			if(inherits(trash2, "try-error")){
			  # warn that the attempt to get ordination coordinates for ordering failed.
			  warning("Attempt to access ordination coordinates for feature/species/taxa/OTU ordering failed.\n",
			          "Using default feature/species/taxa/OTU ordering.")
			}
			if(!is.null(specDF)){
			  # If the score accession seemed to work, go ahead and replace sample.order
			  taxa.order = taxa_names(physeq)[order(phyloseq:::RadialTheta(specDF))]
			}
		}
	}
  
  # Now that index orders are determined, check/assign edges of axes, if specified
  if( !is.null(first.sample) ){
    sample.order = phyloseq:::chunkReOrder(sample.order, first.sample)
  }
  if( !is.null(first.taxa) ){
    taxa.order = phyloseq:::chunkReOrder(taxa.order, first.taxa)
  }

	# melt physeq with the standard user-accessible data melting function
	# for creating plot-ready data.frames, psmelt.
	# This is also used inside some of the other plot_* functions.
	adf = psmelt(physeq)	
	# Coerce the main axis variables to character. 
	# Will reset them to factor if re-ordering is needed.
	adf$OTU = as(adf$OTU, "character")
	adf$Sample = as(adf$Sample, "character")
	if( !is.null(sample.order) ){
		# If sample-order is available, coerce to factor with special level-order
		adf$Sample = factor(adf$Sample, levels=sample.order)
	} else {
		# Make sure it is a factor, but with default order/levels
		adf$Sample = factor(adf$Sample)
	}
	if( !is.null(taxa.order) ){
		# If OTU-order is available, coerce to factor with special level-order
		adf$OTU = factor(adf$OTU, levels=taxa.order)
	} else {
		# Make sure it is a factor, but with default order/levels
		adf$OTU = factor(adf$OTU)
	}

	## Now the plotting part
	# Initialize p.
	p = ggplot(adf, aes(x = Sample, y = OTU, fill=Abundance)) + 
    geom_raster()

	# # Don't render labels if more than max.label
	# Samples
	if( nsamples(physeq) <= max.label ){
		# Add resize layer for samples if there are fewer than max.label
		p <- p + theme(
			axis.text.x = element_text(
				size=phyloseq:::manytextsize(nsamples(physeq), 4, 30, 12),
				angle=-90, vjust=0.5, hjust=0
			)
		)		
	} else {
	  # Remove the labels from any rendering.
		p = p + scale_x_discrete("Sample", labels="")
	}
	# OTUs
	if( ntaxa(physeq) <= max.label ){
		# Add resize layer for OTUs if there are fewer than max.label
		p <- p + theme(
			axis.text.y = element_text(
				size=phyloseq:::manytextsize(ntaxa(physeq), 4, 30, 12)
			)
		)
	} else {
		# Remove the labels from any rendering.
		p = p + scale_y_discrete("OTU", labels="")
	}
	
	# # Axis Relabeling (Skipped if more than max.label):
	# Re-write sample-labels to some sample variable...
	if( !is.null(sample.label) & nsamples(physeq) <= max.label){
		# Make a sample-named char-vector of the values for sample.label
		labvec = as(get_variable(physeq, sample.label), "character")
		names(labvec) <- sample_names(physeq)
		if( !is.null(sample.order) ){
			# Re-order according to sample.order
			labvec = labvec[sample.order]			
		}
		# Replace any NA (missing) values with "" instead. Avoid recycling labels.
		labvec[is.na(labvec)] <- ""
		# Add the sample.label re-labeling layer
		p = p + scale_x_discrete(sample.label, labels=labvec)
	}
	if( !is.null(taxa.label) & ntaxa(physeq) <= max.label){
		# Make a OTU-named vector of the values for taxa.label
		labvec <- as(tax_table(physeq)[, taxa.label], "character")
		names(labvec) <- taxa_names(physeq)
		if( !is.null(taxa.order) ){		
			# Re-order according to taxa.order
			labvec <- labvec[taxa.order]
		}
		# Replace any NA (missing) values with "" instead. Avoid recycling labels.
		labvec[is.na(labvec)] <- ""		
		# Add the taxa.label re-labeling layer
		p = p + scale_y_discrete(taxa.label, labels=labvec)
	}
	
	# Color scale transformations
	if( !is.null(trans) ){
		p = p + scale_fill_gradient(low=low, high=high, trans=trans, na.value=na.value)
	} else {
		p = p + scale_fill_gradient(low=low, high=high, na.value=na.value)	
	}
	
	# Optionally add a title to the plot
	if( !is.null(title) ){
		p = p + ggtitle(title)
	}
			
	return(p)
}
