% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genetic-load-gw.R
\name{aggregate_multiple_consequence_variants}
\alias{aggregate_multiple_consequence_variants}
\title{Process variants with multiple consequences}
\usage{
aggregate_multiple_consequence_variants(consequence_results)
}
\arguments{
\item{consequence_results}{A tibble the output consequence summary of the
function \code{genetic_load_FUNSEQ_and_Consequence_summary}.}
}
\value{
A tibble with same format as input
}
\description{
Some consequences have the following format:
\code{"splice_region_variant&intron_variant&non_coding_transcript_variant"}
This function select the variants that has the highest IMPACT, the IMPACT
rank is HIGH > MODERATE > LOW > MOFIFIER. If there are ties, then the first
variant is selected.
}
