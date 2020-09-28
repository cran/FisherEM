n = 100

test_that("test Chang1983 simulation", {
  simu = simu_bfem(n, which = "Chang1983")
  expect_equal(nrow(simu$Y), n)
  expect_equal(ncol(simu$Y), 15)
  expect_equal(length(simu$cls), n)
  expect_equal(length(unique(simu$cls)), 2)
})

test_that("test section4.2 simulation", {
p = 25
noise = 3
simu = simu_bfem(n, which = "section4.2", p = p, noise = noise)
  expect_equal(nrow(simu$Y), n)
  expect_equal(ncol(simu$Y), p)
  expect_equal(length(simu$cls), n)
  expect_equal(length(unique(simu$cls)), 3)
})
# 
# 
# test_that("test section4.3 simulation", {
#   snr = 10
#   simu = simu_bfem(n, which = "section4.3", snr = 10)
#   expect_equal(nrow(simu$Y), n)
#   expect_equal(ncol(simu$Y), p)
#   expect_equal(length(simu$cls), n)
#   expect_equal(length(unique(simu$cls)), 3)
# })
