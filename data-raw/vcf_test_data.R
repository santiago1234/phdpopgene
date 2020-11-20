## code to prepare test_vcf data
## this dataset corresponds to the first 300
## variants present in chromosome 2
## the purpose of this data is for testing and development

library(VariantAnnotation)
library(tidyverse)


chrn <- "22"
file_vcf_geno <- "~/Dropbox/proyectos/201015-RotationProject/data/201018-PrepareInputDataForRfmix/vcf/query.chr22.vcf.gz"

chunk_size <- 300
tab_geno <- VcfFile(file_vcf_geno, yieldSize = chunk_size)
open(tab_geno)
test_vcf <- readVcf(tab_geno, "hg38")
close(tab_geno)

# load corresponding variant annotation for this region in chr 22

file_vcf_varannotaion <- "~/Dropbox/proyectos/201015-RotationProject/data/201115-VariantAnnotation1k/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.annotation.vcf.gz"

min_range <- test_vcf %>%
  rowRanges() %>%
  ranges() %>%
  start() %>%
  min()
max_range <- test_vcf %>%
  rowRanges() %>%
  ranges() %>%
  end() %>%
  max()
var_range <- IRanges(start = min_range, end = max_range)
rng <- GRanges(seqnames = chrn, ranges = var_range)

tab_annotation <- TabixFile(file_vcf_varannotaion)
svp <- ScanVcfParam(which = rng)

test_anno_vcf <- readVcf(tab_annotation, "hg38", param = svp)

usethis::use_data(test_vcf, overwrite = TRUE)
usethis::use_data(test_anno_vcf, overwrite = TRUE)
