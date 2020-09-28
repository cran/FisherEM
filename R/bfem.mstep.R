bfem.mstep = function(Y, U, tau, Varmeank, Varcovk, model) {
  # 12 different submodels: [DkBk] ... [AkjBk]
  # inspired from fem.mtep.R with the addition of variational correciton in BFEM  
  
  # Initialization
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  K = ncol(tau)
  d = ncol(U)
  U = as.matrix(U)
  
  mu   = matrix(NA,K,d)
  m   = matrix(NA,K,p)
  prop = rep(c(NA),1,K)
  D = array(0,c(K,d,d))
  b = rep(NA,K)
  
  # Projection
  X = as.matrix(Y) %*% as.matrix(U)
  
  # Estimation
  test = 0
  bk = 0
  for (k in 1:K){
    nk  = sum(tau[,k])
    # if (nk ==0) stop("some classes become empty\n",call.=FALSE)
    # if (nk ==0) return(NULL)
    # Prior Probability
    prop[k] = nk / n
    # Mean in the latent space
    mu[k,]  = Varmeank[,k]
    
    m[k,]  = mu[k,] %*% t(U)
    YY  = as.matrix(Y - t(m[k,]*matrix(1,p,n)))
    if (nk < 1) denk = 1 else denk = nk
    YYk = tau[,k] * matrix(1,n,p) * YY
    
    # Estimation of Delta k amongst 12 submodels
    if (model %in% c('DkBk', 'DkB', 'DBk', 'DB')) {
      D[k,,] = crossprod(YYk%*%U, YY%*%U) / denk
    }
    else if (model %in% c('AkjBk', 'AkjB', 'AjBk', 'AjB')) {
      D[k,,] = diag(diag(crossprod(YYk %*%U, YY%*%U) / denk), d)
    }
    else if (model %in% c('AkBk', 'AkB', 'ABk', 'AB')) {
      D[k,,] = diag(rep(Trace(crossprod(YYk %*%U, YY%*%U) / denk)/d,d), d)
    }
    
    
    # compute b[k]
    if (model %in% c('DkB','DB','AkjB', 'AjB', 'AkB','AB')) {
      bk =  bk + (1 / n) * (sum(YYk*YY) - Trace(D[k,,]))
    } else {
      b[k] = (sum(YYk*YY) - Trace(D[k,,])) / (denk * (p-d))
    }
    
    # add variational correction
    if (model %in% c('DkBk', 'DkB', 'DBk', 'DB')) {
      D[k,,] = D[k,,] + Varcovk[,,k]
    }
    else if (model %in% c('AkjBk', 'AkjB', 'AjBk', 'AjB')) {
      D[k,,] = D[k,,] + diag(diag(as.matrix(Varcovk[,,k])), d)
    } else {
      D[k,,] = D[k,,] + diag(Trace(Varcovk[,,k])/d, d)
    }
  }
  
  # Compute Sigma after loop for homoscedastic models
  Sigma = matrix(0, d, d) 
  if (model %in% c('DBk', 'DB', 'AjBk', 'AjB', 'AB', 'ABk')) {
    for (k in 1:K) {
      nk = sum(tau[,k])
      if (nk < 1) denk = 1 else denk = nk
      Sigma = Sigma + (denk / n) * D[k,,]
    }
    for (k in 1:K) {
      D[k,,] = Sigma
    }
  }
  
  
  # compute beta for '*B' models
  if (model %in% c('DkB','DB','AkjB','AkB','AjB','AB')) {
    bk = bk / (p-d)
    b = rep(bk, K)
  }
  
  # avoid numerical
  b[b<=0] = 1e-3
  for (k in 1:K) if (Trace(D[k,,]<1e-3)!=0) test = test+1
  
  
  
  Sk = list()
  Sk$Sigmak = Sk$invSigmak = array(0, dim = c(d,d,K))
  Sk$logdetk = rep(NA, K)
  for (k in 1:K) {
    Sk$Sigmak[,,k] = D[k,,]
    Sk$invSigmak[,,k] = solve(Sk$Sigmak[,,k])
    Sk$logdetk[k] = c(determinant(as.matrix(Sk$Sigmak[,,k]))$modulus)
  }
  Sk$Beta = b
  
  return(list(Sk=Sk, PI=prop))
}
