% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr-verbs.R
\name{select-phyloseq}
\alias{select-phyloseq}
\alias{select_tax_table}
\alias{select_tax_table,taxonomyTable-method}
\alias{select_tax_table,phyloseq-method}
\alias{select_sample_data}
\alias{select_sample_data,sample_data-method}
\alias{select_sample_data,phyloseq-method}
\title{Subset columns in the taxonomy table or sample data using their names and types}
\usage{
select_tax_table(x, ...)

\S4method{select_tax_table}{taxonomyTable}(x, ...)

\S4method{select_tax_table}{phyloseq}(x, ...)

select_sample_data(x, ...)

\S4method{select_sample_data}{sample_data}(x, ...)

\S4method{select_sample_data}{phyloseq}(x, ...)
}
\arguments{
\item{x}{A \code{phyloseq}, \code{taxonomyTable}, or \code{sample_data} object}

\item{...}{Expressions passed to \code{dplyr::select()}}
}
\description{
These functions are wrappers around \code{dplyr::select()} that provide
convenient ways to modify \code{tax_table(x)} and \code{sample_data(x)}.
See \code{dplyr::select()} for supported syntax and helpers.
}
\details{
The special column names '.otu' and '.sample' should not be used; see
\code{mutate-phyloseq} for the ability to change taxa and sample names using
these names.
}
\examples{
data(GlobalPatterns)

GlobalPatterns \%>\% rank_names
GlobalPatterns \%>\% sample_variables
ps <- GlobalPatterns \%>\%
  select_tax_table(Phylum, Genus:Species) \%>\%
  select_sample_data(!dplyr::contains("Barcode"))
ps \%>\% rank_names
ps \%>\% sample_variables
}
