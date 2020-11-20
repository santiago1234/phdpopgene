test_that("count_alt_alleles", {
  expect_equal(1, count_alt_alleles("1|0"))
  expect_equal(0, count_alt_alleles("0|0"))
  expect_equal(2, count_alt_alleles("1|1"))
})
