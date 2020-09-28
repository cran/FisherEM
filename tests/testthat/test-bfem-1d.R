set.seed(12345)
n = 100
simu = simu_bfem(n, which = "Chang1983")
Y = simu$Y
p = ncol(Y)
K = 2
d = K - 1

test_that("bfem kmeans init", {
  res.bfem = bfem(Y, K, model="DB", init = 'kmeans', nstart = 2, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("bfem random init", {
  res.bfem = bfem(Y, K, model="DB", init = 'random', nstart = 2, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("bfem user init", {
  Tinit = t(rmultinom(n, 1, prob = rep(1/2,2)))
  res.bfem = bfem(Y, K, model="DB", init = 'user', Tinit = Tinit, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("test all models", {
  res.bfem = bfem(Y, K, model="all", init = 'kmeans', nstart = 1, mc.cores = 2)
  expect_equal(res.bfem$K, K)
  expect_equal(res.bfem$d, d)
  expect_equal(dim(res.bfem$P), c(n, K))
  expect_equal(dim(res.bfem$var_param$Varmeank), c(d, K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(d, d, K))
  expect_equal(dim(res.bfem$U), c(p, d))
  expect_equal(dim(res.bfem$proj), c(n, d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})

test_that("test 3 models, K grid", {
  res.bfem = bfem(Y, K = 2:6, model=c('DkBk', 'AkjBk', 'AB'), init = 'kmeans', nstart = 1, 
                  maxit.em = 10, eps.em = 1e-3, maxit.ve = 3, mc.cores = 2)
  expect_equal(dim(res.bfem$var_param$Varmeank), c(res.bfem$d, res.bfem$K))
  expect_equal(dim(res.bfem$var_param$Varcovk), c(res.bfem$d, res.bfem$d, res.bfem$K))
  expect_equal(dim(res.bfem$U), c(p, res.bfem$d))
  expect_equal(dim(res.bfem$proj), c(n, res.bfem$d))
  expect_equal(length(res.bfem$elbos), res.bfem$n_ite+1)
})
