\name{fem}
\alias{fem}

\title{
The Fisher-EM algorithm
}
\description{
The Fisher-EM algorithm is a subspace clustering method for high-dimensional data. It is based on the Gaussian Mixture Model and on the idea that the data lives in a common and low dimensional subspace. An EM-like algorithm estimates both the discriminative subspace and the parameters of the mixture model.
}
\usage{
fem(Y,K=2:6,model='AkjBk',method='gs',crit='icl',maxit=50,eps=1e-4,init='kmeans',
                nstart=5,Tinit=c(),kernel='',disp=FALSE,mc.cores=(detectCores()-1),
                subset=NULL)
}
\arguments{
  \item{Y}{ 
The data matrix. Categorical variables and missing values are not allowed.
}
  \item{K}{
An integer vector specifying the numbers of mixture components (clusters) among which the model selection criterion will choose the most appropriate number of groups. Default is 2:6.
}
 \item{model}{A vector of discriminative latent mixture (DLM) models to fit. There are 12 different models: "DkBk", "DkB", "DBk", "DB", "AkjBk", "AkjB", "AkBk", "AkBk", "AjBk", "AjB", "ABk", "AB".  The option "all" executes the Fisher-EM algorithm on the 12 DLM models and select the best model according to the maximum value obtained by model selection criterion. 
} 
 \item{method}{
The method used for the fitting of the projection matrix associated to the discriminative subspace. Three methods are available: 'gs' (Gram-Schmidt, the original proposition), 'svd' (based on SVD, fastest approach, it should be preferred on large data sets) and 'reg' (the Fisher criterion is rewritten as a regression problem). The 'gs' method is the default method since it is the most efficient one on most data sets.
}
   \item{crit}{The model selection criterion to use for selecting the most appropriate model for the data. There are 3 possibilities: "bic", "aic" or "icl". Default is "icl".
}
  \item{maxit}{
The maximum number of iterations before the stop of the Fisher-EM algorithm. 
}
  \item{eps}{
The threshold value for the likelihood differences to stop the Fisher-EM algorithm.
}
  \item{init}{
The initialization method for the Fisher-EM algorithm. There are 4 options: "random" for a randomized initialization, "kmeans" for an initialization by the kmeans algorithm, "hclust" for hierarchical clustering initialization or "user" for a specific initialization through the parameter "Tinit". Default is "kmeans". Notice that for "kmeans" and "random", several initializations are asked and the initialization associated with the highest likelihood is kept (see "nstart").
}
 \item{nstart}{The number of restart if the initialization is "kmeans" or "random". In such a case, the initialization associated with the highest likelihood is kept.
}
  \item{Tinit}{
A n x K matrix which contains posterior probabilities for initializing the algorithm (each line corresponds to an individual).
}
  \item{kernel}{It enables to deal with the n < p problem. By default, no kernel (" ") is used. But the user has the choice between 3 options for the kernel: "linear", "sigmoid" or "rbf".
}
  \item{disp}{
If true, some messages are printed during the clustering. Default is false.
}
  \item{mc.cores}{The number of CPUs to use to fit in parallel the different models (only for non-Windows platforms). Default is the number of available cores minus 1.}

\item{subset}{A positive integer defining the size of the subsample, default is NULL. In case of large data sets, it might be useful to fit a FisherEM model on a subsample of the data, and then use this model to predict cluster assignments for the whole data set. Notice that in, such a case, likelihood values and model selection criteria are computed for the subsample and not the whole data set.}
}

\value{
A list is returned: 
  \item{K}{The number of groups.}
  \item{cls}{the group membership of each individual estimated by the Fisher-EM algorithm.}
  \item{P}{the posterior probabilities of each individual for each group.}
  \item{U}{The loading matrix which determines the orientation of the discriminative subspace.}
  \item{mean}{The estimated mean in the subspace.}
  \item{my}{The estimated mean in the observation space.}
  \item{prop}{The estimated mixture proportion.}
  \item{D}{The covariance matrices in the subspace.}
  \item{aic}{The value of the Akaike information criterion.}
  \item{bic}{The value of the Bayesian information criterion.}
  \item{icl}{The value of the integrated completed likelihood criterion.}
  \item{loglik}{The log-likelihood values computed at each iteration of the FEM algorithm.}
  \item{ll}{the log-likelihood value obtained at the last iteration of the FEM algorithm.}
  \item{method}{The method used.}
  \item{call}{The call of the function.}
  \item{plot}{Some information to pass to the plot.fem function.}
  \item{crit}{The model selction criterion used.}
}
\references{Charles Bouveyron and Camille Brunet (2012), Simultaneous model-based clustering and visualization in the Fisher discriminative subspace, Statistics and Computing, 22(1), 301-324 <doi:10.1007/s11222-011-9249-9>.

Charles Bouveyron and Camille Brunet (2014), "Discriminative variable selection for clustering with the sparse Fisher-EM algorithm", Computational Statistics, vol. 29(3-4), pp. 489-513 <10.1007/s00180-013-0433-6>.
}

\author{
Charles Bouveyron, Camille Brunet & Nicolas Jouvin.
}

\seealso{sfem, plot.fem, fem.ari, summary.fem}

\examples{
data(iris)
res = fem(iris[,-5],K=3,model='AkBk',method='gs')
res
plot(res)
fem.ari(res,as.numeric(iris$Species))
table(iris$Species,res$cls)

\donttest{
# Fit several models and numbers of groups (use by default on non-Windows
# platforms the parallel computing).
res = fem(iris[,-5],K=2:6,model='all',method='gs', mc.cores=2)
res
plot(res)
fem.ari(res,as.numeric(iris$Species))
table(iris$Species,res$cls)
}
}
