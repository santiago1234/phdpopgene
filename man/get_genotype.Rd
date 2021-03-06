% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utlis-genetic-load.R
\name{get_genotype}
\alias{get_genotype}
\title{Extract genotypes from VCF in a tidy frame}
\usage{
get_genotype(vcf_genotypes)
}
\arguments{
\item{vcf_genotypes}{CollapsedVCF, genotypes vcf}
}
\value{
A data frame
\describe{
\item{varid}{Variant ID}
\item{sample1}{Genotype for individual 1}
...
\item{sampleN}{Genotype for individual N}
}
}
\description{
Extract genotypes from VCF in a tidy frame
}
