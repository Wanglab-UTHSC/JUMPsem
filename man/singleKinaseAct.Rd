% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/singleKinaseAct.R
\name{singleKinaseAct}
\alias{singleKinaseAct}
\title{singleKinaseAct}
\usage{
singleKinaseAct(kinase, input, adj, cor.off, kmo.off, mdsite)
}
\arguments{
\item{kinase}{Single kinase name.}

\item{input}{Normalized and transformed phosphoproteomics  data.}

\item{adj}{Adjacency matrix.}

\item{cor.off}{Set up correlation cutoff value 0-1 to remove high collinear variables. Default is 0.95.}

\item{kmo.off}{Set up KMO cutoff value 0-1. Default is 0.}

\item{mdsite}{Logical. Mapping with precise phosphorylated modification sites or not. Default is TRUE. description}
}
\description{
Get single kinase activity.
}
