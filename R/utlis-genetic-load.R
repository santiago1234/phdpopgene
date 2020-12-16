#' Counts the number of alternative alleles
#'
#' @param gtype The genotype, example, 0|0, 1|0, etc.
#'
#' @return number of alternative alleles
#'
count_alt_alleles <- function(gtype) {
  stopifnot(stringr::str_detect(gtype, "\\|"))
  allele_1 <- stringr::str_sub(gtype, 1, 1) %>% as.integer()
  allele_2 <- stringr::str_sub(gtype, 3, 3) %>% as.integer()

  c(allele_1, allele_2) %>%
    (function(x) x > 0)() %>%
    sum()
}


#' Counts the number of alternative alleles
#'
#' A funseq score greater than 1.5 is considered deleterious. An indicator function
#' that returns true if the variant is deleterious multiplied by the number
#' of alternative alleles in the genotype.
#'
#' @param alt_alleles int, how many alternative alles are in the genotype?,
#' for example 0|0 is 0, 1|0 is 1,
#' @param funseqscore the funseq score
#' @param cut_off_val, The value used, FUNSEQ score, to classify a variant as deleterious
#' @return a number 0, 1, or 2.
funseq_load <- function(alt_alleles, funseqscore, cut_off_val=1.5) {
  (funseqscore > cut_off_val) * alt_alleles
}



#' Extract genotypes from VCF in a tidy frame
#'
#' @param vcf_genotypes CollapsedVCF, genotypes vcf
#'
#' @return A data frame
#' \describe{
#'   \item{varid}{Variant ID}
#'   \item{sample1}{Genotype for individual 1}
#'   ...
#'   \item{sampleN}{Genotype for individual N}
#' }
#'
get_genotype <- function(vcf_genotypes) {
  VariantAnnotation::geno(vcf_genotypes)$GT %>%
    tidyr::as_tibble(rownames = "varid")
}



#' Merge genotypes with annotation
#'
#' Merges the genotypes with the annotation (FUNSEQ and Consequences)
#' and puts the table in long format. This function also removes
#' all sites, for each individual-sample, that are homozygous (0|0) since
#' these sites don't contribute to the load.
#'
#' @param genotipos the tidy genotypes, output of \code{\link{get_genotype}}
#' @param h_samples vector of samples, a subset of the samples
#' @param annotation table with annotation, the output of \code{\link{get_var_annotation}}
#' @return A data frame, with the following columns
#' \describe{
#'   \item{varid}{Variant ID}
#'   \item{Consequence}{Variant consequence}
#'   \item{FUNSEQ}{funseq score}
#'   \item{individual}{individual id, example: HG00551}
#'   \item{genotype}{the genotype of the individual, example: 0|1}
#' }
merge_genotypes_with_annotation <- function(genotipos, annotation, h_samples) {
  all(h_samples %in% colnames(genotipos)) %>% stopifnot()

  genotipos[, c("varid", h_samples)] %>%
    dplyr::inner_join(annotation, by = "varid") %>%
    tidyr::pivot_longer(
      cols = -c(.data$varid, .data$FUNSEQ, .data$Consequence),
      names_to = "individual", values_to = "genotype"
    ) %>%
    dplyr::filter(.data$genotype != "0|0") %>%
    dplyr::mutate(n_alternative = purrr::map_int(.data$genotype, .f = count_alt_alleles))

}


#' Summarize the Consequence
#'
#' For each individual computes a summary of how many variants types (e. g.
#' intergenic_variant, stop_lost, etc.) carries.
#'
#' Consider a particular variant and suppose the consequence is stop_gained.
#' If the individual is 0|0 it wont carrie the variant, if it is 0|1 or 0|1, it
#' will carry the variance 1 and 2 for 1|1
#' @param geno_anno_mlong Tibble with genotypes and annotation, the output
#' of \code{\link{merge_genotypes_with_annotation}}
#'
#' @return A data frame, with the following columns
#' \describe{
#'   \item{individual}{individual id, example: HG00551}
#'   \item{Consequence}{Variant consequence}
#'   \item{n_variants}{How many variants of type Consequence the}
#' }
consequence_summary <- function(geno_anno_mlong) {
  geno_anno_mlong %>%
    dplyr::group_by(.data$individual, .data$Consequence) %>%
    dplyr::summarise(n_variants = sum(.data$n_alternative)) %>%
    dplyr::ungroup()
}


#' How many variants with (FUNSEQ > cut_off_val) do each individual carry?
#'
#' For each individual computes how many variants with FUNSEQ > cut_off_val
#' carries.
#'
#' If the individual is 0|0 we add nothing to funseq_load
#' If the individual is 1|0 or 0|1 we add 1 to funseq_load
#' If the individual is 1|1 or 0|1 we add 2 to funseq_load
#' @inheritParams consequence_summary
#' @param cut_off_val, The value used, FUNSEQ score, to classify a variant as deleterious
#'
#' @return A data frame, with the following columns
#' \describe{
#'   \item{individual}{individual id, example: HG00551}
#'   \item{funseq_load}{How many variants have a funseqscore > cut_off_val}
#' }
funq_load_summary <- function(geno_anno_mlong, cut_off_val) {
  geno_anno_mlong %>%
    dplyr::mutate(
      fl = purrr::map2_int(.x = .data$n_alternative, .y = .data$FUNSEQ, .f = funseq_load, cut_off_val = cut_off_val)
    ) %>%
    dplyr::group_by(.data$individual) %>%
    dplyr::summarise(funseq_load = sum(.data$fl, na.rm = TRUE))
}



#' Computes FUNSEQ load and Consequence summary.
#'
#' The funseq load for an individual is defined as: TODO
#'
#' NOTE: If an individual has non variant sites (i. e. 0|0) it wont be
#' present in the results table.
#'
#' @inheritParams get_var_annotation
#' @param cores int, the number of cores used for parallel computation. The
#' computation is parallelized over the individuals (samples).
#' @param cut_off_val, The value used, FUNSEQ score, to classify a variant as deleterious
#' @return A list with two elements: FS_l (funseq load) and CSQ_s
#' @export
#' @example
#' res <- genetic_load_FUNSEQ_and_Consequence_summary(test_vcf, test_anno_vcf,22, 1)
#' cqs <- res$CSQ_s
#' funseq <- res$FS_l
genetic_load_FUNSEQ_and_Consequence_summary <- function(vcf_genotypes, vcf_annotation, chr, cores = 6, cut_off_val = 1.5) {

  # sanity check
  # check given chromosome matches range in vcf files
  rng <- SummarizedExperiment::rowRanges(vcf_genotypes)
  chr_found <- GenomicRanges::seqnames(rng) %>%
    as.character() %>%
    unique()

  if (!chr_found == chr) {
    stop("Supplied chromosome does not match chromosome in VCF file")
  }



  genotipos <- get_genotype(vcf_genotypes)
  annotation <- get_var_annotation(vcf_genotypes, vcf_annotation, chr)

  message("Computing load ...")
  # helper functions for parallel computation
  fs_load <- function(h_ss) {
    m_l <- merge_genotypes_with_annotation(genotipos, annotation, h_ss)
    funq_load_summary(m_l, cut_off_val = cut_off_val)

  }
  csq_summary <- function(h_ss) {
    m_l <- merge_genotypes_with_annotation(genotipos, annotation, h_ss)
    consequence_summary(m_l)

  }


  # split sample to perform parallel computation
  h_samples <- colnames(genotipos)
  h_samples <- h_samples[2:length(h_samples)]

  # split into n groups of size 10 each
  h_samples <- split(h_samples, ceiling(seq_along(h_samples) / 10))

  f_l <- parallel::mclapply(h_samples, fs_load, mc.cores = cores) %>%
    purrr::reduce(dplyr::bind_rows)

  c_s <- parallel::mclapply(h_samples, csq_summary, mc.cores = cores) %>%
    purrr::reduce(dplyr::bind_rows)


  list(
    FS_l = f_l,
    CSQ_s = c_s
  )


}
