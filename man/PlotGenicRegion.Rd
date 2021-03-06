% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting_functions.R
\name{PlotGenicRegion}
\alias{PlotGenicRegion}
\title{Plot genic region}
\usage{
PlotGenicRegion(genic.dt, type = c("stacked", "bar"), filename)
}
\arguments{
\item{genic.dt}{hybrid data.table}

\item{type}{type of plot, either stacked or bar}

\item{filename}{filename for plot with extension}
}
\value{
plot
}
\description{
Plot the hybrid region for all arms for genic reads and facet by intra/inter
}
