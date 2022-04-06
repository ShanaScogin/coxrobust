

test_that("coxr works", {
  # Create a simple test data set using the attached function gen_data
  set.seed(123)
  a <- gen_data(200, c(1, 0.1, 2), cont = 0.05, p.censor = 0.30)
  result <- coxr(Surv(time, status) ~ X1 + X2 + X3, data = a , trunc = 0.9)
  coefs <- unlist(result[][1])
  names(coefs) <- c("1", "2", "3")
  check_against <- c(0.752665504, 0.208303005, 1.498068463)
  names(check_against) <- c("1", "2", "3")
  expect_equal(coefs, check_against)
})


test_that("coxr works", {
  # Create a simple test data set using the attached function gen_data
  set.seed(123)
  a <- gen_data(200, c(1, 0.1, 2), cont = 0.05, p.censor = 0.30)
  result <- coxr(Surv(time, status) ~ X1 + X2 + X3, data = a , trunc = 0.9)
  coefs <- unlist(result[][1])
  names(coefs) <- c("1", "2", "3")
  check_against <- c(0.752665504, 0.208303005, 1.498068463)
  names(check_against) <- c("1", "2", "3")
  expect_equal(coefs, check_against)
})
