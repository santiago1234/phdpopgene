# a set of function to parse vcf annotation from 1k genomes


#' Extract FUNSEQ score from vcf annotation info
#'
#' extract the funseq scores in a tidy data frame from 1k variant annottation.
#' The FunSeq score is an evaluation of the deleterious effect of noncoding variants.
#' FunSeq score > 1.5 is suggested to be a cutoff for being deleterious variants,
#' and used here in some analyses
#' \url{https://www.internationalgenome.org/faq/where-can-i-get-consequence-annotations-1000-genome-variants/}
#' @param vcf_annotation CollapsedVCF, object read with \code{VariantAnnotation::readVcf}
#'
#' @return A data frame
#' \describe{
#'   \item{varid}{Variant ID}
#'   \item{FUNSEQ}{funseq score}
#' }
#' @export
#'
#' @examples
#' get_FUNSEQ(test_anno_vcf)
get_FUNSEQ <- function(vcf_annotation) {
  VariantAnnotation::info(vcf_annotation) %>%
    tidyr::as_tibble(rownames = "varid") %>%
    dplyr::select(.data$varid, .data$FUNSEQ) %>%
    dplyr::mutate(
      FUNSEQ = purrr::map_dbl(.data$FUNSEQ, function(x) x[[1]]) # extract 1st element
    )
}


#' Extract Calculated variant consequence
#'
#' get the Calculated variant consequence
#' \url{https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html}
#' @param vcf_annotation CollapsedVCF, object read with \code{VariantAnnotation::readVcf}
#'
#' @return A data frame
#' \describe{
#'   \item{varid}{Variant ID}
#'   \item{Consequence}{Calculated variant consequence}
#' }
#' @export
#'
#' @examples
#' get_CSQ(test_anno_vcf)
get_CSQ <- function(vcf_annotation) {

  # from the VCF file header we have:
  # CSQ,Number=.,Type=String,Description="Consequence type as predicted by VEP WITH -PICK_ALLELE parameter. Format: Allele|Gene|Feature|Feature_type|Consequence|...
  csq_names <-
    "Allele|Gene|Feature|Feature_type|Consequence|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|SIFT|PolyPhen|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE" %>%
    stringr::str_split(pattern = "\\|") %>%
    base::unlist()

  VariantAnnotation::info(vcf_annotation) %>%
    tidyr::as_tibble(rownames = "varid") %>%
    dplyr::select(.data$varid, .data$CSQ) %>%
    dplyr::mutate(CSQ = purrr::map_chr(.data$CSQ, function(x) x[[1]])) %>%
    tidyr::separate(.data$CSQ, into = csq_names, sep = "\\|") %>%
    dplyr::select(.data$varid, .data$Consequence)
}



#' Add rs ID to INFO annotation
#'
#' the vcf file with variant annotation sometimes does not
#' include the rs IDs and returns another id when I extract the info
#' this new id has the format: chrn:start_RF/ALT for some variants
#' @param vcf_genotypes CollapsedVCF, genotypes vcf
#' @param info_data this can be the output of \code{\link{get_CSQ}} or \code{\link{get_FUNSEQ}}
#' @param chrn int or chr, chromosome name
#' @return A data frame
#' \describe{
#'   \item{varid}{Variant ID}
#'   \item{range_id}{Id for variant in range chrn:start_RF/ALT}
#'   ...
#' }
#' @export
#'
#' @examples
#' info_funseq <- get_FUNSEQ(test_anno_vcf)
#' add_rsids_to_funseq(test_vcf, info_funseq, 22)
add_rsids_to_funseq <- function(vcf_genotypes, info_data, chrn) {

  # which variant ids are nat an rsID and start with chrn?
  ids_chr <- info_data %>%
    dplyr::filter(stringr::str_detect(.data$varid, paste0("^", chrn))) %>%
    dplyr::rename(range_id = .data$varid)

  # which are rs ids
  ids_rs <- info_data %>%
    dplyr::filter(
      !stringr::str_detect(.data$varid, paste0("^", chrn))
    )

  # generate a table mapping chromosome position id to rs-ids

  ranges <- SummarizedExperiment::rowRanges(vcf_genotypes)
  # add a column with the var-id

  ranges$varid <- names(ranges)
  map_of_ids <-
    unique(ranges) %>%
    as.data.frame() %>%
    tidyr::as_tibble() %>%
    dplyr::mutate(
      ALT = purrr::map_chr(.data$ALT, function(x) as.character(x[[1]])), # character because the ALT is a Biostring
      range_id = paste(.data$seqnames, ":",
        .data$start, "_",
        .data$REF, "/",
        .data$ALT,
        sep = ""
      )
    ) %>%
    dplyr::select(.data$varid, .data$range_id)

  dplyr::inner_join(map_of_ids, ids_chr, by = "range_id") %>%
    dplyr::bind_rows(
      dplyr::inner_join(map_of_ids, ids_rs, by = "varid")
    )
}


#' Extract variant annotation
#'
#' This function extracts the variant annotation FUNSEQ and Consequence from
#' the annotation vcf data
#'
#' @inheritParams add_rsids_to_funseq
#' @param vcf_annotation CollapsedVCF -> \code{VariantAnnotation::readVcf}. The annotaion for
#' the variants in \code{vcf_genotypes}
#' @return A data frame
#' \describe{
#'   \item{varid}{Variant ID}
#'   \item{range_id}{Id for variant in range chrn:start_RF/ALT}
#'   \item{Consequence}{Calculated variant consequence}
#'   \item{FUNSEQ}{FUNSEQ score}
#' }
#' @export
#'
#' @examples
#' get_var_annotation(test_vcf, test_anno_vcf, 22)
get_var_annotation <- function(vcf_genotypes, vcf_annotation, chrn) {

  funseq_scores <- get_FUNSEQ(vcf_annotation)
  funseq_scores <- add_rsids_to_funseq(vcf_genotypes, funseq_scores, chrn) %>%
    dplyr::select(-.data$range_id)

  csq <- get_CSQ(vcf_annotation)
  csq <- add_rsids_to_funseq(vcf_genotypes, csq, chrn) %>%
    dplyr::select(-.data$range_id)

  dplyr::inner_join(csq, funseq_scores, by = "varid")
}
