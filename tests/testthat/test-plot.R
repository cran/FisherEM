library(ggplot2)
set.seed(12345)
n= 200
simu = simu_bfem(n, which = "Chang1983")
Y = simu$Y
p = ncol(Y)
K = 2
d= K - 1
res.bfem = bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2)

test_that("Plot bound", {
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "elbo")))
})

test_that("Plot crit", {
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "criterion", crit = 'bic')))
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "criterion", crit = 'icl')))
  expect_true(is.ggplot(plot.bfem(res.bfem, type = "criterion", crit = 'aic')))
})

# d == 1 non implemented yet
# test_that("Plot subspace", {
#   expect_true(is.ggplot(plot.bfem(res.bfem, type = "subspace")))
# })