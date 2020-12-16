## quiets concerns of R CMD check re: the .'s that appear in pipelines
utils::globalVariables(c(".", ".data"))
