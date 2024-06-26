% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singleEnzymeAff.R
\name{singleEnzymeAff}
\alias{singleEnzymeAff}
\title{singleEnzymeAff}
\usage{
singleEnzymeAff(enzyme, input, adj, cor.off, kmo.off, mdsite)
}
\arguments{
\item{enzyme}{Single enzyme name.}

\item{input}{Normalized and transformed phospho-proteomics / ubiqutin-proteomics data.}

\item{adj}{Adjacency matrix.}

\item{cor.off}{Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.}

\item{kmo.off}{Set up KMO cutoff value 0-1. Default is 0.}

\item{mdsite}{Logical. Mapping with precise phosphorylated modification sites or not. Default is TRUE.}
}
\description{
Get single enzyme affinity.
}
