#' Experimental setting of the chapter BFEM
#'
#' @param n Number of observations
#' @param which Type of simulation, either:
#' \itemize{
#' \item "Chang1983" - Simulate the dataset of Chang's (1983) paper : a mixture of 2 Gaussian with in dimension p=15. 
#' \item "section4.2" - Experimental setting of Section 4.2: DLM model in dimension p with d=2 and K=3, with noisy dimensions. 
#' \item "section4.3" - Experimental setting of Section 4.3: Same as `"section4.2"` except the noise is expressed in term of signal-to-noise ration (decibels).
#' }
#' @param ... Additional param controlling the simulation
#' \itemize{
#' \item p - The desired observed space dimension, the latent dimension is kept fixed to d=2 and noisy Gaussian dimensions are added (useless for `"Chang1983"`) 
#' \item noise (for `"section4.2"` only) - Variance of the noise 
#' \item snr (for `"section4.3"` only) - Signal-to-noise ratio (in decibels) representing the ratio of signal and noise variances in logarithmic scale. The greater snr, the smaller noise variance.
#' }
#' 
#' @return A list with slots
#' \itemize{
#' \item Y - The simulated data.
#' \item cls - The true clustering.
#' }
#' 
#' @export
#' @examples
#' n = 300
#' 
#' # Chang's 1983 setting
#' simu = simu_bfem(n = n, which = "Chang1983")
#' 
#' # Section 4.2 setting
#' p = 25
#' noise = 1
#' simu = simu_bfem(n, which = "section4.2", p = p, noise = noise)
#' 
#' # Section4.3 setting
#' snr = 3 # noise variance is 2 times smaller than that of the signal.
#' simu = simu_bfem(n, which = "section4.3", snr = 10)
simu_bfem <- function(n, which = "Chang1983", ...) {
  simu = switch(which,
                "Chang1983" = simu.Chang1983(n=n),
                "section4.2" = simu.section4.2(n=n, ...),
                "section4.3" = simu.section4.3(n=n, ...))
  return(simu)
}

# Chang 1983 dataset
simu.Chang1983 <- function(n = 300) {
  p = 15
  
  cte = matrix(rep(0.95 - 0.05* (1:p), n), nrow = n, byrow = T)
  Z = stats::rbinom(n=n, size=1, prob = 0.2)
  
  mu = rep(0, p)
  Sigma = diag(1, p)
  f = c(rep(-0.9, 8), rep(0.5, 7))
  for(i in 1:p) {
    for (j in 1:p) {
      if (j!= i) Sigma[i,j] = -0.13 * f[i] * f[j] 
    }
  }
  X = MASS::mvrnorm(n=n, mu = mu, Sigma = Sigma)
  
  Y = 0.5 * cte + Z * cte + X
  return(list(Y=Y, cls = factor(Z)))
}

# High-dimensional setting of Section 4.2
simu.section4.2 <- function(n, p = 50, noise = 1) {
  d = 2
  K = 3
  S = array(0,c(d,d,K))
  S[,,2] = S[,,3] = S[,,1] = matrix(c(1.5,0.75, 0.75, 0.45))
  
  cst = 2.5
  mu = matrix(0, d, K)
  mu[,1] = c(0, 0)
  mu[,2] = cst * c(1, 0)
  mu[,3] = cst * c(2, 0)

  n1 = round(4*n/10)
  n2 = round(3*n/10)
  n3 = round(3*n/10)
  
  ind = sample(1:n)
  X = rbind(mvrnorm(n1,mu[,1],S[,,1]),mvrnorm(n2,mu[,2],S[,,2]),mvrnorm(n3,mu[,3],S[,,3]))
  X = X[ind,]
  A = mvrnorm(p, mu = rep(0,p),Sigma = diag(25, p, p))
  W = qr.Q(qr(A))
  B = mvrnorm(n,rep(0, p-d), diag(noise, p-d, p-d))
  y = cbind(X,B)
  Y = y %*% t(W)
  
  cls = c(rep(1,n1),rep(2,n2),rep(3,n3))[ind]
  
  return(list(Y = Y, W = W, X = X, cls = cls))  
}

# Signal-to-noise ratio of section 4.3
simu.section4.3 <- function(n, p = 150, snr = 3) {
  # inverse formula: snr = 10*log10(1.95 / beta)
  beta = 1.95 / (exp(snr * log(10)/10))
  simu.section4.2(n=n, p = p, noise = beta)
}