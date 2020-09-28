bfem.main <- function(Y,K,init,nstart,control_bfem,Tinit,model,kernel,method,lambda){

  em.tol = control_bfem$em$tol
  ve.tol = control_bfem$var$tol
  ve.max.iter = control_bfem$var$max.iter
  em.max.iter = control_bfem$em$max.iter
  emp.bayes = control_bfem$emp.bayes
  
  # Initialization
  colNames = colnames(Y)
  Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  d = min((K-1),(p-1))
  
  # Compute S
  m = colMeans(Y)
  XX = as.matrix(Y - t(m*t(matrix(1,n,p))))
  S = t(XX) %*% XX /n

  
  # ========== Initialization of tau with a frequentist EM-step ==============
  if (init=='user'){ tau = Tinit}
  else if (init=='hclust'){
    tau   = matrix(0,n,K)
    ind = cutree(hclust(dist(Y),method='ward.D2'),K)
    for (i in 1:n){ tau[i,ind[i]] = 1 }
  }
  else if (init=='kmeans' || init=='random'){
    Ltmp = rep(NA,nstart); tau_list = list()
    for (i in 1:nstart){
      if (init=='random'){tau_list[[i]] = t(rmultinom(n,1,c(rep(1/K,K))))}
      else{
        tau = matrix(0,n,K)
        ind = kmeans(Y,K,nstart=10)$cluster
        for (i in 1:n){ tau[i,ind[i]] = 1 }
        tau_list[[i]] = tau
      }
      U = switch(method,
                 'svd'= fstep.fisher(XX,tau_list[[i]],S,kernel),
                 'gs'= fstep.GramSc(XX,tau_list[[i]],S,kernel),
                 'reg'  = fstep.qiao(Y,tau_list[[i]],kernel))
      prms      = fem.mstep(Y,U,tau_list[[i]],model=model,method=method)
      res.estep = fem.estep(prms,Y,U)
      Ltmp[i]   = res.estep$loglik
    }
    tau = tau_list[[which.max(Ltmp)]]
  }
  Tinit = tau # stock to return
  logtau = log(tau)
  
  
  # ================ Set initial quantities fo BFEM ==========================
  # initial F-step
  U = switch(method,'svd'= fstep.fisher(XX,tau,S,kernel),
             'gs'= fstep.GramSc(XX,tau,S,kernel),
             'reg'  = fstep.qiao(Y,tau,kernel))
  
  # Initialize nu
  if(emp.bayes) {
    X = Y %*% U
    nu = colMeans(X)
  } else {
    nu = rep(0, d)
  }
  
  # Initial M-step: use frequentist estimates of Sk and PI
  prms.fem = fem.mstep(Y,U,tau,model=model,method=method)
  prms = transform_param_fem_to_bfem(prms.fem)
  
  # Inital VE-step: 1 iteration of fixed point for var parameter of q(\mu_k)
  Varcovk = updateVarcovk(prms, lambda, tau)
  Varmeank = updateVarmeank(Y, U, prms, nu, tau, Varcovk)
  
  elbo = bfem.elbo(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk)
  
  # ================ Main loop BFEM ===============================
  conv = c(Inf)
  Lobs = rep(-Inf , 1, em.max.iter + 1)
  Lobs[1] = elbo
  Linf_new  = Lobs[1]
  for (i in 1:em.max.iter) {
    # if (verbose > 0) cat('BFEM iteration : ', i,'\n')
    
    # F-step
    U = switch(method,
               'svd'= fstep.fisher(XX,tau,S,kernel),
               'gs'= fstep.GramSc(XX,tau,S,kernel),
               'reg'  = fstep.qiao(Y,tau,kernel))
    
    # Variational E-step (don't change tau yet, only q(\mu_k))
    ve_step = bfem.vestep(Y, U, prms, nu, lambda, logtau, 
                          Varmeank=NULL, Varcovk=NULL, ve.tol, ve.max.iter)
    Varmeank = ve_step$Varmeank
    Varcovk = ve_step$Varcovk
    
    # Mstep
    prms = bfem.mstep(Y, U, tau, Varmeank, Varcovk, model)

    if(is.null(prms$Sk)) return(list(k=K, Tinit=Tinit, loglik=-Inf, 
                                icl = -Inf, bic = -Inf, aic = -Inf, 
                                model = model, K=K, cls = max.col(tau), P = tau,
                                error=1))
    
    
    # Update of tau and logtau after M-step
    logtau = updatelogTau(Y, U, prms, Varmeank, Varcovk) 
    tau = exp(logtau)
    
    # Hyper-params 
    if (emp.bayes) {
      nu = updateNu(Varmeank)
      lambda = updateLambda(Varmeank, Varcovk, nu)
    }
    
    elbo = bfem.elbo(Y, U, prms, nu, lambda, logtau, Varmeank, Varcovk)
    # Stop criterion
    if (is.na(elbo)) {
      warning('Likelihood becomes NA (probably emptied cluster).')
      break
    }
    Lobs[i+1] = elbo
    
    if (i >= 3) {
      acc = (Lobs[i+1] - Lobs[i]) / (Lobs[i] - Lobs[i-1])
      Linf_old = Linf_new
      Linf_new <- try( Lobs[i] + 1/(1 - acc) * (Lobs[i+1] - Lobs[i]))
      if (is.na(Linf_new)){warnings('acc=', acc);break}
      if (abs(Linf_new - Linf_old) < em.tol) {break}
    } 
  }
  Lobs = Lobs[Lobs != -Inf]
  
  
  # Returning the results
  cls = apply(tau, 1, which.max) 
  param = list(PI = prms$PI, Sigmak = prms$Sk$Sigmak, Beta = prms$Sk$Beta)
  var_param = list(logtau = logtau,  Varmeank = Varmeank, Varcovk = Varcovk)
  hypers = list(lambda = lambda, nu = nu)
  proj = Y %*% U
  rownames(U) = colNames
  final_elbo = utils::tail(Lobs, 1)
  crit = bfem.criteria(Y, final_elbo, U, prms, tau, model, hypers, var_param, method)
  return(list(K = K, Tinit = Tinit, d = d, cls = cls, P = tau,
              elbos = Lobs, loglik=final_elbo, n_ite = i,
              U = U, param = param, var_param = var_param, proj = proj,
              model = model, aic=crit$aic, bic=crit$bic, icl=crit$icl,
              comp=crit$comp,
              hypers = hypers,
              call = call)
         )
}

transform_param_fem_to_bfem <- function(prms.fem) {
  # New code use different structure for parameters
  # This a helper function to transform prms.fem
  # into the structure used in bfem.main()
  
  K = prms.fem$K
  p = prms.fem$p
  d = min(p-1, K-1)
  
  # new structure
  prms = list(Sk = list(), PI = NULL)
  prms$PI = prms.fem$prop
  prms$Sk = list(Sigmak = array(0, dim = c(d,d,K)),
                 invSigmak = array(0, dim = c(d,d,K)),
                 Beta = prms.fem$b,
                 logdetk = rep(NA, K))
  for (k in 1:K) {
    prms$Sk$Sigmak[,,k] = prms.fem$D[k,,]
    prms$Sk$invSigmak[,,k] = solve(prms$Sk$Sigmak[,,k])
    prms$Sk$logdetk[k] = c(determinant(as.matrix(prms$Sk$Sigmak[,,k]))$modulus)
  }
  return(prms)
}
