% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/please.R
\name{summarizeAndPlotMutations}
\alias{summarizeAndPlotMutations}
\title{Summarize and Plot Mutation Types with Context}
\usage{
summarizeAndPlotMutations(
  vcf_file,
  ref_genome,
  context_length = 3,
  plot = TRUE,
  top_n = 20
)
}
\arguments{
\item{vcf_file}{File path to the VCF file or a CollapsedVCF object.}

\item{ref_genome}{Reference genome object.}

\item{context_length}{Integer specifying the length of the genomic context (default: 3).}

\item{plot}{Logical indicating whether to plot the mutation counts (default: TRUE).}

\item{top_n}{Integer specifying the number of most common mutation types to display in the plot (default: 20).}
}
\value{
A data frame containing the counts of each mutation type with context.
}
\description{
Summarizes the mutation types extracted from a VCF file and plots
the frequency of each mutation type categorized by reference and mutated nucleotides.
It displays a refined plot showing only the most common mutation types for legibility,
while returning a complete frequency table.
}
\examples{
summarizeAndPlotMutations(vcf, Hsapiens,3)

}
