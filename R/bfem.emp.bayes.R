
updateNu = function(Varmeank) {
  return(rowMeans(Varmeank))
}

updateLambda = function(Varmeank, Varcovk, nu) {
  d = dim(Varcovk)[1]
  K = dim(Varcovk)[3]
  
  lambda = 0
  for (k in 1:K) {
    lambda = lambda + Trace(Varcovk[,,k]) + sum(Varmeank[,k]^2) 
    - 2 * t(nu) %*% Varmeank[,k] 
    + sum(nu^2)
  }
  return(lambda / (d * K))
}
