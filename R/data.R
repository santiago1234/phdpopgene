#' A subset of vcf file for chromosome 22 (1000 genomes)
#'
#' VCF file used for testing function and package development
#' @format A VCF file with 300 variants and a subset of samples from America
#' @source \url{https://github.com/santiago1234/phdpopgene/blob/main/data-raw/vcf_test_data.R}
"test_vcf"


#' Variant annotation
#'
#' the variant annotation for \code{\link{test_vcf}}
#' @format A VCF file with 300 variants and its corresponding annotation
#' @source \url{https://github.com/santiago1234/phdpopgene/blob/main/data-raw/vcf_test_data.R}
"test_anno_vcf"


#' Ensembl Variation - Calculated variant consequences
#'
#' @format A data frame with 36 rows and 3 variables:
#' \describe{
#'   \item{Consequence}{Variant Consequence}
#'   \item{Display_term}{Variant Consequence in a nice format}
#'   \item{IMPACT}{Variant Impact, HIGH means deleterious}
#'   ...
#' }
#' @source \url{https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html}
"variant_consequences"
