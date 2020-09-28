bfem.vestep <- function(Y, U, prms, nu, lambda, logtau, 
                        Varmeank = NULL, Varcovk = NULL, 
                        ve.tol = 1e-4, ve.max.iter = 3) {
  
  tau = exp(logtau)
  ve.count = 0 
  if (is.null(Varmeank) | is.null(Varcovk)) {
    Varcovk = updateVarcovk(prms, lambda, tau)
    Varmeank = updateVarmeank(Y, U, prms, nu, tau, Varcovk)
  }
  
  elbo = bfem.elbo(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk)
  ve.elbos = c(elbo)
  conv = c(elbo, Inf)
  test = F
  
  while (abs(diff(conv)/conv[1]) > ve.tol & ve.count < ve.max.iter) {
    
    logtau = updatelogTau(Y, U, prms, Varmeank, Varcovk)
    tau = exp(logtau)
    
    Varcovk = updateVarcovk(prms, lambda, tau)
    Varmeank = updateVarmeank(Y, U, prms, nu, tau, Varcovk) 
    
    elbo = bfem.elbo(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk)
    if (is.na(elbo)) {test = T; break}
    ve.elbos = c(ve.elbos, elbo)
    conv = c(elbo, conv[1])
    ve.count = ve.count + 1
  }
  if (sum(diff(ve.elbos) >= -1e-10) != length(ve.elbos) - 1) {
    message('Elbo decrease : ', diff(ve.elbos))
    # warning(paste0('The elbo is decreasing in VE-step'))
  }
  
  return(list(logtau = logtau, Varmeank = Varmeank, Varcovk = Varcovk, 
              ve.elbos = ve.elbos, n_iter = ve.count, test = test))
}

updatelogTau <- function(Y, U, prms, Varmeank, Varcovk) {
  # Initialization
  Sk = prms$Sk
  PI = prms$PI
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  d = ncol(U)
  K = ncol(Varmeank)
  prop = PI
  D = array(0, dim = c(K, d, d))
  for (k in 1:K)  D[k,,] = Sk$Sigmak[,,k]
  b = Sk$Beta
  QQ = matrix(NA,n,K)
  T = matrix(NA,n,K)
  
  # Compute posterior probabilities
  for (k in 1:K){
    bk = b[k]
    mY = t(U %*% Varmeank[,k])
    YY = Y-t(matrix(rep(mY,n),p,n)) 
    projYY = YY %*% U %*% t(U)
    if (d==1){
      QQ[,k] =  (1/D[k,1,1]) * rowSums(projYY^2) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(D[k,1,1]) - 2*log(prop[k]) + p*log(2*pi)
      # add variational correction
      QQ[,k] = QQ[,k] + (1/D[k,1,1]) * Varcovk[1,1,k]
      
    } else{
      tmp = eigen(D[k,(1:d),(1:d)])
      A = projYY %*% U %*% tmp$vect %*% diag(sqrt(1/tmp$val))
      QQ[,k] = rowSums(A^2) + 1/bk*rowSums((YY - projYY)^2) + (p-d)*log(bk) + log(det(D[k,(1:d),(1:d)])) - 2*log(prop[k]) + p*log(2*pi)
      # add variational correction
      QQ[,k] = QQ[,k] + sum(Varcovk[,,k] * Sk$invSigmak[,,k])
      
    }
  }
  # Compute posterior probabilities
  for (k in 1:K) {T[,k] = 1 / rowSums(exp((QQ[,k]*matrix(1,n,K)-QQ)/2))}
  
  # Return the results
  return(log(T))
}



updateVarcovk = function(prms, lambda, tau) {
  # tau is a matrix with dim (nxK)
  
  Sk = prms$Sk
  nk = colSums(tau)
  K = ncol(tau)
  d = dim(Sk$invSigmak)[1]
  res = array(0, dim = c(d,d, K))
  for (k in 1:K) {
    res[,,k] = solve(diag(1/lambda, nrow = d) +  nk[k] * Sk$invSigmak[,,k])
  }
  return(res)
}

updateVarmeank = function(Y, U, prms, nu, tau, Varcovk) {
  # Varcovk is an array with dim (d x d x K)
  # invSigmak is an array with dim (p x p x K)
  
  Sk = prms$Sk
  K = ncol(tau)
  nk = colSums(tau)
  d = ncol(U)
  
  res = matrix(0, nrow = d, ncol = K)
  for (k in 1:K) {
    sumYpondk = colSums(tau[,k] * Y)
    res[,k] = nu + Varcovk[,,k] %*% Sk$invSigmak[,,k] %*% (t(U) %*% sumYpondk - nk[k] * nu )
    
  }
  
  return(res)
}

