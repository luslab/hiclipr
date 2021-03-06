% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hybrid_functions.R
\name{FindNonIslandHybrids}
\alias{FindNonIslandHybrids}
\title{Find non-island hybrid reads}
\usage{
FindNonIslandHybrids(reads, islands)
}
\arguments{
\item{reads}{hybrid data.table of reads}

\item{islands}{hybrid data.table of islands}
}
\value{
hybrid data.table
}
\description{
Find hybrid reads that do not overlap an island
}
