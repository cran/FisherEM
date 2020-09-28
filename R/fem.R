fem <- function(Y,K=2:6,model='AkjBk',method='gs',crit='icl',maxit=50,eps=1e-4,init='kmeans',
                nstart=5,Tinit=c(),kernel='',disp=FALSE,mc.cores=(detectCores()-1),subset=NULL){
  call = match.call()
  MOD = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB','AjB','AjBk', 'ABk', 'AB','all')
  MET = c('svd','reg','gs')
  KER = c('','sigmoid','linear','rbf')
  CRI = c('bic','aic','icl','sh')
  INIT = c('user','random','kmeans','hclust')
  if (sum(is.na(Y))>0) stop("NA values are not allowed.\n",call.=FALSE)
  if (any(!model%in%MOD)) stop("Invalid model name.\n",call.=FALSE)
  if (any(!method%in%MET)) stop("Invalid method name.\n",call.=FALSE)  
  if (any(!kernel%in%KER)) stop("Invalid kernel name.\n",call.=FALSE)
  if (any(!crit%in%CRI)) stop("Invalid criterion name.\n",call.=FALSE)
  if (any(!init%in%INIT)) stop("Invalid initialization name.\n",call.=FALSE)
  if (init=='hclust' & nrow(Y)>5000) stop('Too much data for this initialization.\n',call.=FALSE)
  # if (K>=ncol(Y)) stop("K must be strictly less than the number of variables",call.=FALSE)
  if (nrow(Y)<=ncol(Y) & method=='gs') stop("n<<p case: use method REG or SVD.\n",call.=FALSE)
  if (length(model)==1) if (model=='all') model = MOD[MOD!='all']
  if (sum(apply(Y,2,var) == 0) > 0) stop("Some variables have zero variance. Please remove them and try again.\n",call.=FALSE)
  
  if (!is.null(subset)) {
    Yfull = Y
    sel =  sample(nrow(Y),subset)
    Y = Y[sel,]
    if (init=='user') Tinit = Tinit[sel,]
  }
  
  # Run FEM depending on Windows or not (allows parallel computing)
  if (Sys.info()[['sysname']] == 'Windows' | mc.cores == 1){
    prms = expand.grid(model=model,K=K)
    RES = list()
    for (i in 1:nrow(prms)){
      RES[[i]] = fem.main(Y=Y,K=prms$K[i],model=prms$model[i],init=init,nstart=nstart,maxit=maxit,
                          eps=eps,Tinit=Tinit,kernel=kernel,method=method)
    }
  }
  else {
    prms = expand.grid(model=model,K=K)
    MoreArgs = list(Y=Y,init=init,nstart=nstart,maxit=maxit,eps=eps,Tinit=Tinit,kernel=kernel,method=method)
    RES = do.call(mcmapply, c(list(FUN="fem.main",MoreArgs=MoreArgs,mc.cores=mc.cores,
                                   mc.silent=TRUE,mc.preschedule=FALSE),prms))
  }
    
  #Post-treatment of results
  if (is.matrix(RES)){ # Parallization without errors (output is a matrix)
    bic = unlist(apply(RES,2,function(x){if (is.list(x)){x$bic} else NA})) 
    aic = unlist(apply(RES,2,function(x){if (is.list(x)){x$bic} else NA}))
    icl = unlist(apply(RES,2,function(x){if (is.list(x)){x$bic} else NA}))
    comp = unlist(apply(RES,2,function(x){if (is.list(x)){x$comp} else NA}))
    loglik = unlist(apply(RES,2,function(x){if (is.list(x)){x$loglik} else NA}))
    if (crit=='bic'){ id_max = which.max(bic); crit_max = bic[id_max]}
    if (crit=='aic'){ id_max = which.max(aic); crit_max = aic[id_max]}
    if (crit=='icl'){ id_max = which.max(icl); crit_max = icl[id_max]}
  } 
  else{ # Parallization with errors (output is a list)
    bic = unlist(lapply(RES,function(x){if(is.list(x)){x$bic} else{-Inf}}))
    aic = unlist(lapply(RES,function(x){if(is.list(x)){x$aic} else{-Inf}}))
    icl = unlist(lapply(RES,function(x){if(is.list(x)){x$icl} else{-Inf}}))
    comp = unlist(lapply(RES,function(x){if(is.list(x)){x$comp} else{-Inf}}))
    loglik = unlist(lapply(RES,function(x){if(is.list(x)){x$loglik} else{-Inf}}))
    if (crit=='bic'){ id_max = which.max(bic); crit_max = bic[id_max]}
    if (crit=='aic'){ id_max = which.max(aic); crit_max = aic[id_max]}
    if (crit=='icl'){ id_max = which.max(icl); crit_max = icl[id_max]}
  }
  nm = length(model)
  allCriteria = data.frame(model=prms$model,K=prms$K,comp=comp,loglik=loglik,bic=bic,aic=aic,icl=icl)
  if (is.matrix(RES)) {res = RES[,id_max]}  else res = RES[[id_max]]
  res$aic = res$bic = res$icl = NULL
  res$model =  as.character(res$model)
  res$allCriteria = allCriteria
  res$crit = crit
  res$critValue = unlist(crit_max)
  res$allResults = RES
  res$call = call
  
  if (!is.null(subset)) {
    browser()
    prms = list(mean=res$mean,my=res$my,K=res$K,prop=res$prop,D=res$D,b=res$b)
    resFull = fem.estep(prms,Yfull,res$U)
    res$P = resFull$T
    res$cls  = max.col(resFull$T)
  }
    
  # Display and return results
  if (disp) cat('The selected model is',res$model,'with K =',res$K,'(',crit,'=',res$critValue,')\n')
  class(res)='fem'
  res
}