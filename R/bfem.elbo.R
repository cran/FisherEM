bfem.elbo <- function(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk) {
  
  # Set quantities 
  Sk = prms$Sk
  PI =  prms$PI
  
  p = ncol(Y)
  K = length(PI)
  d = ncol(U)
  
  tau = exp(logtau)
  tau = tau / rowSums(tau)  
  nk = colSums(tau)
  
  # H(q(z))
  E1 = -sum(tau * logtau, na.rm = T)
  #  E_z[log p(z)]
  E3 = sum(nk * log(PI))
  E2 = E4 = E5 = 0
  for (k in 1:K) {
    
    # H[q(mu)]
    E2 = E2 + 0.5 * (d * (log(2*pi) + 1) + c(determinant(as.matrix(Varcovk[,,k]))$modulus))
    
    # E_mu[log p(mu)]
    
    E4 = E4 - 0.5 * (d * log(2*pi) 
                     + d * log(lambda) 
                     +  (1/lambda) * (Trace(Varcovk[,,k]) + sum(Varmeank[,k]^2) 
                                      - 2 * t(nu) %*% Varmeank[,k] 
                                      + sum(nu^2)
                     )
    )
    
    
    
    # E_{z, mu}[log p(Y|Z,mu)]
    Ypondk = tau[,k] * Y
    sumYpondk = colSums(Ypondk)
    projk = U %*% Varmeank[,k]
    if (nk[k] < 1e-16) {
      E5 = NaN
      warning('Cluster ', k, ' is emptied.')
      break
    } else {
      Ck = (1 / nk[k]) * (t(Ypondk) %*% Y -
                            sumYpondk %*% t(projk) - 
                            projk %*% t(sumYpondk)
      ) +
        U %*% (Varcovk[,,k] + Varmeank[,k] %*% t(Varmeank[,k])) %*% t(U)
      tUCkU = t(U) %*% Ck %*% U
      E5 = E5 - 0.5 * nk[k] * (p * log(2*pi) + Sk$logdetk[k] +
                                 (p - d) * log(Sk$Beta[k]) +
                                 Trace(Sk$invSigmak[,,k] %*% tUCkU) +
                                 (1/Sk$Beta[k]) * (Trace(Ck) - Trace(tUCkU))
      )
      
    }
  }
  bound = c(E1 + E2 + E3 + E4 + E5)
  return(bound)
  
}