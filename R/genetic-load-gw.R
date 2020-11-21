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


#' TODO
#'
#' FUNSEQ load? TODO
#'
#' CONSEQUENCE summary? TODO
#'
#' @param vcf_file
#' @param annotation_file
#' @param chr
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
gl_funseq_consequence_summary <- function(vcf_file, annotation_file, chr, cores = 4) {

  chunk_size <- 1000 # load 1000 SNPs in each iteration

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
      cores = cores
    )

    # add results to tibble
    funseq_load_result <- aggregate_res_fs(funseq_load_result, results$FS_l)
    consequence_s_result <- aggregate_res_cqs(consequence_s_result, results$CSQ_s)

  }

  close(tab_geno)

  message("<- DONE ->")

  list(
    FS_load = funseq_load_result,
    CQS_summary = consequence_s_result
  )

}

# resultados <- gl_funseq_consequence_summary(
#   vcf_file = "~/Dropbox/proyectos/201015-RotationProject/data/201018-PrepareInputDataForRfmix/vcf/query.chr22.vcf.gz",
#   annotation_file = "~/Dropbox/proyectos/201015-RotationProject/data/201115-VariantAnnotation1k/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz",
#   chr = 22
# )






