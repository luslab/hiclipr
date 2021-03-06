% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/barcodes_mapping_functions.R
\name{AlignReadsBowtie}
\alias{AlignReadsBowtie}
\title{Align and merge reads using Bowtie (as in original paper, except no RC for step 3)}
\usage{
AlignReadsBowtie(fq.file, bam.file, processors)
}
\arguments{
\item{fq.file}{FASTQ file (unzipped)}

\item{bam.file}{output BAM filename}

\item{processors}{Number of processors to use}
}
\value{
SAM file of alignments
}
\description{
Align reads and merge as directed
}
