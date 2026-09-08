test_that("GATES inference agrees with evalITR through the public D interface", {
  set.seed(86)
  n <- 600
  D <- rep(c(0, 1), n / 2)
  Y <- stats::rnorm(n)
  score <- stats::runif(n)
  tau <- matrix(score, n, 3)
  ind <- rep(1:3, each = 200)
  expect_equal(GATE(D, score, Y, 3),
               evalITR::GATE(D, score, Y, 3))
  expect_equal(GATEcv(D = D, tau = tau, Y = Y, ind = ind, ngates = 3),
               evalITR::GATEcv(T = D, tau = tau, Y = Y, ind = ind, ngates = 3))
  expect_equal(hetcv.test(D, tau, Y, ind, 3),
               evalITR::hetcv.test(D, tau, Y, ind, 3))
  expect_equal(het.test(D, score, Y, 3),
               evalITR::het.test(D, score, Y, 3))
  set.seed(14)
  a <- consistcv.test(D, tau, Y, ind, 3, nsim = 99)
  set.seed(14)
  b <- evalITR::consistcv.test(D, tau, Y, ind, 3, nsim = 99)
  expect_equal(a, b)
  set.seed(14)
  a <- consist.test(D, score, Y, 3, nsim = 99)
  set.seed(14)
  b <- evalITR::consist.test(D, score, Y, 3, nsim = 99)
  expect_equal(a, b)
})
