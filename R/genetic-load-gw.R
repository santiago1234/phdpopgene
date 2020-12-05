# stupid test cases
# chromosomes should mach



get_annotation_currentRange <- function(vcf_yield, tab_annotation, chr) {

  # what is the range of vcf_yield ? (Genomic Positions)

  min_range <- vcf_yield %>%
    SummarizedExperiment::rowRanges() %>%
    GenomicRanges::ranges() %>%
    GenomicRanges::start() %>%
    min()

  max_range <- vcf_yield %>%
    SummarizedExperiment::rowRanges() %>%
    GenomicRanges::ranges() %>%
    GenomicRanges::start() %>%
    max()

  var_range <- IRanges::IRanges(start = min_range, end = max_range)
  rng <- GenomicRanges::GRanges(seqnames = chr, ranges = var_range)
  svp <- VariantAnnotation::ScanVcfParam(which = rng)
  VariantAnnotation::readVcf(tab_annotation, "hg38", param = svp)

}


aggregate_res_fs <- function(previous_res, current_res) {
  # add current results to previous

  dplyr::bind_rows(previous_res, current_res) %>%
    dplyr::group_by(.data$individual) %>%
    dplyr::summarise(funseq_load = sum(.data$funseq_load)) %>%
    dplyr::ungroup()
}

aggregate_res_cqs <- function(previous_res, current_res) {
  # add current results to previous

  dplyr::bind_rows(previous_res, current_res) %>%
    dplyr::group_by(.data$individual, .data$Consequence) %>%
    dplyr::summarise(n_variants = sum(.data$n_variants)) %>%
    dplyr::ungroup()
}


#' Process variants with multiple consequences
#'
#' Some consequences have the following format:
#' \code{"splice_region_variant&intron_variant&non_coding_transcript_variant"}
#' This function select the variants that has the highest IMPACT, the IMPACT
#' rank is HIGH > MODERATE > LOW > MOFIFIER. If there are ties, then the first
#' variant is selected.
#'
#' @param consequence_results A tibble the output consequence summary of the
#' function \code{genetic_load_FUNSEQ_and_Consequence_summary}.
#'
#' @return A tibble with same format as input
#'
#' @examples
aggregate_multiple_consequence_variants <- function(consequence_results) {

  # helper functions
  get_var <- function(impact) {
    dplyr::filter(variant_consequences, .data$IMPACT == impact) %>%
      dplyr::pull(.data$Consequence)
  }

  variant_rank <- list(
    HIGH = get_var("HIGH"),
    MODERATE = get_var("MODERATE"),
    LOW = get_var("LOW"),
    MODIFIER = get_var("MODIFIER")
  )

  highest_consequence_var <- function(consecuencias) {
    # consecuencias is a character vector of consequences
    extractor <- function(var_list, from) {
      x_matching <- var_list[var_list %in% from]
      x_matching[1]
    }

    # the order below gives the priority
    if (any(consecuencias %in% variant_rank$HIGH)) {
      return(extractor(consecuencias, variant_rank$HIGH))
    }

    if (any(consecuencias %in% variant_rank$MODERATE)) {
      return(extractor(consecuencias, variant_rank$MODERATE))
    }

    if (any(consecuencias %in% variant_rank$LOW)) {
      return(extractor(consecuencias, variant_rank$LOW))
    }

    if (any(consecuencias %in% variant_rank$MODIFIER)) {
      return(extractor(consecuencias, variant_rank$MODIFIER))
    }
  }

  worst_var <- function(multivar) {
    # multivar is a string that contains multiple variants separated by
    # the character &

    multivar %>%
      stringr::str_split("&") %>%
      unlist() %>%
      highest_consequence_var()
  }

  multiple_consequences <- dplyr::filter(consequence_results, stringr::str_detect(.data$Consequence, "&"))
  unique_consequence <- dplyr::filter(consequence_results, !stringr::str_detect(.data$Consequence, "&"))


  # get worst consequence and aggregate to unique consequence

  multiple_consequences %>%
    dplyr::mutate(
      Consequence = purrr::map_chr(.data$Consequence, worst_var)
    ) %>%
    dplyr::bind_rows(unique_consequence) %>%
    dplyr::group_by(.data$individual, .data$Consequence) %>%
    dplyr::summarise(n_variants = sum(.data$n_variants)) %>%
    dplyr::ungroup()



}


#' TODO
#'
#' FUNSEQ load? TODO
#'
#' CONSEQUENCE summary? TODO
#'
#' @param vcf_file
#' @param annotation_file
#' @param chr
#' @param cores integer, number of cores used for parallel processing
#' @prama cut_off_val, The value used, FUNSEQ score, to classify a variant as deleterious
#'
#' @return
#' @export
#'
#' @examples
gl_funseq_consequence_summary <- function(vcf_file, annotation_file, chr, cores = 4, cut_off_val = 1.5) {

  chunk_size <- 10000 # load 10000 SNPs in each iteration

  tab_geno <- VariantAnnotation::VcfFile(vcf_file, yieldSize = chunk_size)
  tab_annotation <- Rsamtools::TabixFile(annotation_file)

  open(tab_geno)


  current_range <- 0

  # create list to save intermediate results
  funseq_load_result <- tibble::tibble()
  consequence_s_result <- tibble::tibble()

  while (nrow(vcf_yield <- VariantAnnotation::readVcf(tab_geno, "hg38"))) {

    message("processing variants:", current_range, "-", current_range + chunk_size, "\n")
    current_range <- current_range + chunk_size

    # get the annotation for vcf_yield
    annotaion_yield <- get_annotation_currentRange(vcf_yield, tab_annotation, chr)

    # analyze the load in current variants ------------------------------------

    results <- genetic_load_FUNSEQ_and_Consequence_summary(
      vcf_genotypes = vcf_yield,
      vcf_annotation = annotaion_yield,
      chr = chr,
      cores = cores,
      cut_off_val = cut_off_val
    )

    # add results to tibble
    funseq_load_result <- aggregate_res_fs(funseq_load_result, results$FS_l)
    consequence_s_result <- aggregate_res_cqs(consequence_s_result, results$CSQ_s)

  }

  close(tab_geno)

  # additional processing for consequence summary
  # agregate multivars and add column with consequence
  consequence_s_result <-
    aggregate_multiple_consequence_variants(consequence_s_result) %>%
    dplyr::inner_join(variant_consequences, by = "Consequence")

  message("<- DONE ->")

  list(
    FS_load = funseq_load_result,
    CQS_summary = consequence_s_result
  )

}
