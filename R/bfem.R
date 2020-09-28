
#' The Bayesian Fisher-EM algorithm.
#'
#' The Bayesian Fisher-EM algorithm is built on a Bayesian formulation of the
#' model used in the \code{\link{fem}}. It is a subspace clustering method for
#' high-dimensional data. It is based on a Gaussian Mixture Model and on the
#' idea that the data lives in a common and low dimensional subspace. A VEM-like
#' algorithm estimates both the discriminative subspace and the parameters of
#' the mixture model.
#'
#' @param Y The data matrix. 
#' Categorical variables and missing values are not
#'   allowed.
#' @param K An integer vector specifying the numbers of mixture components
#'   (clusters) among which the model selection criterion will choose the most
#'   appropriate number of groups. Default is 2:6.
#' @param model A vector of Bayesian discriminative latent mixture (BDLM) models
#'   to fit. There are 12 different models: "DkBk", "DkB", "DBk", "DB", "AkjBk",
#'   "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB".  The option "all"
#'   executes the Fisher-EM algorithm on the 12 DLM models and select the best
#'   model according to the maximum value obtained by model selection criterion.
#'   Similar to \code{\link{fem}}
#' @param method The method used for the fitting of the projection matrix
#'   associated to the discriminative subspace. Three methods are available:
#'   'gs' (Gram-Schmidt, the original proposition), 'svd' (based on SVD, faster)
#'   and 'reg' (the Fisher criterion is rewritten as a regression problem). The
#'   'gs' method is the default method.
#' @param crit The model selection criterion to use for selecting the most
#'   appropriate model for the data. There are 3 possibilities: "bic", "aic" or
#'   "icl". Default is "icl".
#' @param maxit.em The maximum number of iterations before the stop of the main
#'   EM loop in the BFEM algorithm.
#' @param eps.em The threshold value for the likelihood differences (Aitken's
#'   criterion) to stop the BFEM algorithm.
#' @param maxit.ve The maximum number of iterations before the stop of the
#'   VE-step loop (fixed point algorithm)
#' @param eps.ve The threshold value for the likelihood differences (Aitken's
#'   criterion) to stop the BFEM algorithm.
#' @param lambda The initial value for the variance of the Gaussian prior on the
#'   means in the latent space.
#' @param emp.bayes Should the hyper-parameters (mean and variance) of the prior be updated ? Default to TRUE.
#' @param init The initialization method for the Fisher-EM algorithm. There are
#'   4 options: "random" for a randomized initialization, "kmeans" for an
#'   initialization by the kmeans algorithm, "hclust" for hierarchical
#'   clustering initialization or "user" for a specific initialization through
#'   the parameter "Tinit". Default is "kmeans". Notice that for "kmeans" and
#'   "random", several initializations are asked and the initialization
#'   associated with the highest likelihood is kept (see "nstart").
#' @param nstart The number of restart if the initialization is "kmeans" or
#'   "random". In such a case, the initialization associated with the highest
#'   likelihood is kept.
#' @param Tinit A n x K matrix which contains posterior probabilities for
#'   initializing the algorithm (each line corresponds to an individual).
#' @param kernel It enables to deal with the n < p problem. By default, no
#'   kernel (" ") is used. But the user has the choice between 3 options for the
#'   kernel: "linear", "sigmoid" or "rbf".
#' @param disp If true, some messages are printed during the clustering. Default
#'   is false.
#' @param mc.cores The number of CPUs to use to fit in parallel the different
#'   models (only for non-Windows platforms). Default is the number of available
#'   cores minus 1.
#' @param subset A positive integer defining the size of the subsample, default
#'   is NULL. In case of large data sets, it might be useful to fit a FisherEM
#'   model on a subsample of the data, and then use this model to predict
#'   cluster assignments for the whole data set. Notice that in, such a case,
#'   likelihood values and model selection criteria are computed for the
#'   subsample and not the whole data set.
#'
#' @return A list is returned: \itemize{  
#' \item K - The number of groups.
#' \item cls - the group membership of each individual estimated by the BFEM algorithm 
#' \item Tinit - The initial posterior probalities used to start the algorithm 
#' \item d - the dimension of the discriminative subspace 
#' \item elbos - A vector containing the evolution of the variational lower bound at each iteration 
#' \item loglik - The final value of the variational lower bound 
#' \item n_ite - The number of iteration until convergence of the BFEM algorithm 
#' \item P - the posterior probabilities of each individual for each group 
#' \item U - The loading matrix which determines the orientation of the discriminative subspace 
#' \item param - A list containing the estimated parameters of the model 
#' \itemize{
#' \item PI - The mixture proportions 
#' \item Sigmak - An array containing estimated cluster covariances in the latent space 
#' \item Beta - The noise variance in each cluster 
#' }
#' \item var_param - A list containing the variational distribution parameters
#' \itemize{
#' \item logtau - A n x K matrix containing the logarithm of the multinomial parameters of q(Z) 
#' \item Varmeank - A K x d matrix containing the variational mean 
#' \item Varcovk -  A d x d x K array containing the variational covariance matrices.
#' }
#' \item proj - The projected data on the discriminative subspace. 
#' \item aic - The value of the Akaike information criterion 
#' \item bic - The value of the Bayesian information criterion 
#' \item icl - The value of the integrated completed likelihood criterion 
#' \item method - The method used in the F-step 
#' \item call - The call of the function 
#' \item crit - The model selection criterion used }
#' @export
#' @seealso \code{\link{fem}}
#' @examples
#' # Chang's 1983 setting
#' simu = simu_bfem(300, which = "Chang1983")
#' Y = simu$Y
#' res.bfem = bfem(Y, K = 2:6, model=c('AB'), init = 'kmeans', nstart = 1, 
#'                maxit.em = 10, eps.em = 1e-3, maxit.ve = 3, mc.cores = 2)
#' 
bfem <- function(Y, K=2:6, model='AkjBk', method='gs', crit='icl', maxit.em=100,
                eps.em=1e-6, maxit.ve=3, eps.ve=1e-4, lambda = 1e3, emp.bayes=T,
                init='kmeans', nstart=10, Tinit=c(), kernel='', 
                disp=FALSE, mc.cores=(detectCores()-1), subset=NULL) {
  
  call = match.call()
  control_bfem = list(var = list(tol = eps.ve, max.iter = maxit.ve),
                      em = list(tol = eps.em, max.iter = maxit.em),
                      emp.bayes = emp.bayes)
  MOD = c('DkBk','DkB','DBk','DB','AkjBk','AkjB','AkBk','AkB',
          'AjB','AjBk', 'ABk', 'AB','all')
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
      RES[[i]] = bfem.main(Y=Y,K=prms$K[i],model=prms$model[i],init=init,
                           nstart=nstart, control_bfem=control_bfem,
                           Tinit=Tinit,kernel=kernel,method=method, lambda=lambda)
    }
  }
  else {
    prms = expand.grid(model=model,K=K)
    MoreArgs = list(Y=Y,init=init,nstart=nstart,control_bfem=control_bfem,
                    Tinit=Tinit,kernel=kernel,method=method, lambda=lambda)
    RES = do.call(mcmapply, c(list(FUN="bfem.main",MoreArgs=MoreArgs,mc.cores=mc.cores,
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
  class(res)='bfem'
  res
}