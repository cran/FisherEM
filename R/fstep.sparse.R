fstep.sparse <-
function(X,T,lambda,nbit,l2){

    if (length(lambda)>1) stop('The user needs to enter a single figure comprised between 0 and 1')
# 	require('lars')
# 	require('elasticnet')
	# Initialization
	K = ncol(T)
	p = ncol(X)
	d = min(p-1,(K-1))
	m = matrix(NA,K,p)

	# Compute summary statistics
	Xbar = colMeans(X)
	n = colSums(T)
	for (k in 1:K){ m[k,] = colSums((as.matrix(T[,k]) %*% matrix(1,1,p))* X) / n[k] }
	
	# Matrices Hb and Hw
	Hb =  as.matrix(sqrt(n) * (m - matrix(1,K,1) %*% Xbar))
	Hw = X - t(apply(T,1,'%*%',m))
	
	# Cholesky decomposition of t(Hw) %*% Hw
	if (nrow(X)>p) Rw = chol(t(Hw)%*%Hw) else {
				gamma = 0.5
				Rw = chol(t(Hw)%*%Hw + gamma*diag(p))}

	# LASSO & SVD
	Binit = eigen(ginv(cov(X))%*%(t(Hb)%*%Hb))$vect[,1:d]
	if (is.complex(Binit)) Binit = matrix(Re(Binit),ncol=d,byrow=F)
	if (is.null(dim(Binit))) B = matrix(Binit) else B=Binit
	res.svd = svd(t(ginv(Rw))%*%t(Hb)%*%Hb%*%B)
	A = res.svd$u %*% t(res.svd$v)
	for (i in 1:nbit){
		for (j in 1:d){
			W = rbind(Hb,sqrt(lambda)*Rw)
			y = rbind(Hb %*% ginv(Rw) %*% A[,j],matrix(0,p,1))
# 			res.lasso = lars(W,y,intercept=FALSE)
# 			B[,j] = coef(res.lasso,mode="fraction",s=lambda)
# 			if (n<=p) {l2 = 1}
			res.enet = enet(W,y,lambda=l2,intercept=FALSE)
			B[,j] = predict.enet(res.enet,X,type="coefficients",mode="fraction",s=lambda)$coef
		}
		normtemp = sqrt(apply(B^2, 2, sum))                                                           
		normtemp[normtemp == 0] = 1
		Beta     = t(t(B)/normtemp)
		res.svd  = svd(t(ginv(Rw))%*%t(Hb)%*%Hb%*%Beta)
		A        = res.svd$u %*% t(res.svd$v)
		Beta = svd(Beta)$u
# 		Beta = qr.Q(qr(Beta))
	}
	
	# return the sparse loadings
	as.matrix(Beta)
}

