bfem.criteria = function(Y, final_elbo, U, prms, tau, model, hypers, var_param, method) {
  
  n = nrow(tau)
  K = ncol(tau)
  p = nrow(U)
  d = ncol(U)
  
  nu = hypers$nu
  # set vague prior
  lambda = 1e10
  
  Varmeank = var_param$Varmeank
  Varcovk = var_param$Varcovk
  
  # compute final_elbo with new U and q
  Varcovk = updateVarcovk(prms, lambda, tau)
  Varmeank = updateVarmeank(Y,U, prms, nu, tau = tau, Varcovk = Varcovk)
  prms = bfem.mstep(Y, U, tau, Varmeank, Varcovk, model)
  final_elbo = bfem.elbo(Y, U, prms, nu, lambda, log(tau), Varmeank, Varcovk)
  
  # compute penalty
  if (method=='sparse'){ p = sum(abs(U) > 1e-2) }
  comp = switch(as.character(model),
                'DkBk' = (K-1) + d*(p-(d+1)/2) + K*d*(d+1)/2 + K,
                'DkB'  = (K-1) + d*(p-(d+1)/2) + K*d*(d+1)/2 + 1,
                'DBk'  = (K-1) + d*(p-(d+1)/2) + d*(d+1)/2 + K,
                'DB'   = (K-1) + d*(p-(d+1)/2) + d*(d+1)/2 + 1,
                'AkjBk'= (K-1) + d*(p-(d+1)/2) + K*d + K,
                'AkjB' = (K-1) + d*(p-(d+1)/2) + K*d+1,
                'AkBk' = (K-1) + d*(p-(d+1)/2) + K + K,
                'AkB'  = (K-1) + d*(p-(d+1)/2) + K + 1,
                'AjBk' = (K-1) + d*(p-(d+1)/2) + d + K,
                'AjB'  = (K-1) + d*(p-(d+1)/2) + d + 1,
                'ABk'  = (K-1) + d*(p-(d+1)/2) + 1 + K,
                'AB'   = (K-1) + d*(p-(d+1)/2) + 1 + 1)
  
  if (!is.na(final_elbo)) {
    aic = final_elbo - comp # AIC criterion
    bic = final_elbo - 1/2 * comp * log(n) # BIC criterion
    cl = max.col(tau)
    tau.cl = matrix(0, n, K)
    for (i in 1:n) tau.cl[i, cl[i]] = 1
    icl = final_elbo - 1/2 * comp * log(n) - sum(tau.cl * log(tau+1e-15))
  } else {
    aic = bic = icl = -Inf
  }
  
  list(aic=aic,bic=bic,icl=icl, comp=comp)
}
