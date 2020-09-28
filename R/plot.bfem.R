## For CRAN check ...
## https://stackoverflow.com/questions/9439256/how-can-i-handle-r-cmd-check-no-visible-binding-for-global-variable-notes-when/12429344#12429344
utils::globalVariables(names = c('elbo', 'iteration', 'Cluster',
                                 'x', 'y', '.id', 'sigma', 'V1', 'V2'))


#' Plotting function
#'
#' Utility function to plot the results of the BFEM algorithm. The S3 plot
#' function is a wrapper function over the 3 other functions
#' @param x The results of \code{\link{bfem}}.
#' @param type The plot type: \itemize{ \item "subspace" (default) - Uses
#'   \code{plot_subspace()} to plot the projected data \item "criterion" - Uses
#'   \code{plot_crit()} to plot the criterion value. \item "elbo" - Uses
#'   \code{plot_bound()} to plot the variational lower bound evolution. }
#' @param ... Additional parameter to pass to corxponding functions:
#' @param crit  Used to specify which criterion should be plotted. Possible values are "aic", "bic" and 'icl. The default is the criterion used in the algorithm.
#' @param alpha_levels  A vector giving the desired Gaussian ellipses level set. Default to 0.95.
#' @param plot.dims  The dimension to be plotted. Default to the first two dimensions.
#' @param show.ellipses  Should Gaussian ellipses be plotted. Default to TRUE
#' @param show.uncertainty  Should uncertainty be plotted. A point is considered uncertain if its posterior probability of membership is peaked toward 2 or more clusters. Graphically, it can be displayed with a bigger point size depending on the uncertainty level, bigger points being more uncertain.
#' @param size  The point size.
#' @param cex.uncertainty  The multiplicative factor for the basic point size controlling the size of uncertain points.
#' @return a ggplot2 plot object
#' 
#' @rdname plot.bfem 
#' @export
#' 
#' @examples
#' \donttest{
#' data(iris)
#' Y = iris[,-5]
#' res = bfem(Y, 3, model = 'DB')
#' gg = plot(x=res, type = "subspace")
#' print(gg)
#' }
plot.bfem <- function(x, type = "subspace", ...) {
  gg = switch (type,
               "subspace" = plot_subspace(x, ...),
               "elbo" = plot_bound(x),
               "criterion" = plot_crit(x, ...)
  )
  gg
}


# plot.bfem <- function(x, ...) {
#   mc = match.call()
#   if(is.null(mc$type)) mc$type = "subspace"
#   type = mc$type
#   localPlot_bound <- function(x, ..., type) plot_bound(x, ...)
#   localPlot_crit <- function(x, ..., type) plot_crit(x, ...)
#   localPlot_subspace <- function(x, ..., type) plot_subspace(x, ...)
#   
#   gg = switch (type,
#     "subspace" = localPlot_subspace(x, ...),
#     "elbo" = localPlot_bound(x, ...),
#     "criterion" = localPlot_crit(x, ...)
#   )
#   gg
# }




#--- Plot results from BFEM

bubble <- function(u, cex = c(0.2, 3), alpha = c(0.1, 1)) 
{
  # Uncertainty computing function from Mclust
  # Return size and transparency for points
  u <- as.vector(u)
  cex <- cex[!is.na(cex)]
  alpha <- alpha[!is.na(alpha)]
  u <- (u - min(u))/(max(u) - min(u) + sqrt(.Machine$double.eps))
  n <- length(u)
  r <- sqrt(u/pi)
  r <- (r - min(r, na.rm = TRUE))/
    (max(r, na.rm = TRUE) - min(r, na.rm = TRUE) + sqrt(.Machine$double.eps))
  cex <- r * diff(range(cex)) + min(cex)
  alpha <- u * diff(range(alpha)) + min(alpha)
  return(list(cex = cex, alpha = alpha))
}


#' @describeIn plot.bfem Plot Y projected on the `plot.dims` dimensions of the latent space
# along with gaussian ellipses corxponding to p(\mu_k)
#' @export
plot_subspace <- function(x, alpha_levels = c(0.95), plot.dims = c(1,2), 
                          show.ellipses = T, show.uncertainty = T, 
                          size = 2, cex.uncertainty = 1, ...) {
  
  if (x$d == 1) {
    message('Not implemented for d=1.')
    return(NULL)
  }
  if (sum(plot.dims > x$d) != 0) stop('Plot dimensions must be < d')
  X_est = x$proj
  df = as.data.frame(X_est[,plot.dims])
  df$Cluster = as.factor(x$cl)
  ndim = length(plot.dims)

  if(ndim == 2) {
    gg = ggplot2::ggplot(df, 
                         ggplot2::aes(x = V1, y = V2, col = Cluster, shape = Cluster))
    
    if (show.uncertainty) {
      # inspired from Mclust uncertainty, adapted for ggplot2
      z = x$P
      uncertainty <- 1 - apply(z, 1, max)
      u = (uncertainty - min(uncertainty))/(max(uncertainty) - min(uncertainty)
                                            + sqrt(.Machine$double.eps))
      b <- bubble(u, cex = cex.uncertainty * c(size, size+2), alpha = c(1, 0.7))    
      
      gg = gg + 
        ggplot2::geom_point(size = b$cex,
                   alpha = max(b$alpha) - b$alpha + min(b$alpha))
    } else {
      gg = gg +
        ggplot2::geom_point(size = size)
    }
    
    
    
    # -- code for ellipse
    if (show.ellipses) {
      alpha_levels = sort(alpha_levels)
      names(alpha_levels) <- alpha_levels ## to get id column in xult
      contour_data = NULL
      for (k in 1:x$K) {
        m <- x$var_param$Varmeank[plot.dims,k]
        sigma <- x$param$Sigmak[plot.dims,plot.dims,k]
        contour_data = plyr::ldply(alpha_levels, 
                                   function(level) 
                                     ellipse::ellipse(level = level, x = sigma, centre = m))
        
        contour_data$Cluster = as.factor(k)
        gg = gg +
          ggplot2::geom_path(data=contour_data, 
                             ggplot2::aes(x, y, group = .id, col = Cluster), 
                    linetype='dotdash',
                    alpha = rep(length(alpha_levels):1,each = 100)/1.5 ,
                    show.legend = F)
      }
      
      means = as.data.frame(t(x$var_param$Varmeank))
      means$Cluster = as.factor(1:x$K)
      gg = gg + 
        ggplot2::geom_point(data=means, 
                            ggplot2::aes(x = V1, y = V2),shape = 3, size = 5, show.legend = F)
    }
    
    # Colorblind friendly palette
    gg = gg +
      ggplot2::scale_color_brewer(palette="Set2") +
      ggplot2::xlab(paste0('U', plot.dims[1])) + ggplot2::ylab(paste0('U', plot.dims[2]))
    return(gg)
  } else {
    base::message('ndim > 2 not implemented yet')
    return(NULL)
  }
}

#' @describeIn plot.bfem plot the variational bound evolution 
#' @export
plot_bound = function(x, ...) {
  elbos = x$elbos
  df = data.frame(iteration = 1:length(elbos), elbo = elbos)
  gg = ggplot2::ggplot(df, ggplot2::aes(x=iteration, y=elbo)) +
    ggplot2::geom_point(size=1.5) + 
    ggplot2::geom_line(size=1, alpha = 0.7, linetype='dashed') +
    ggplot2::xlab('Iteration') + ggplot2::ylab('Elbo') +
    ggplot2::theme(text = ggplot2::element_text(size=20)) 
  gg
}

#' @describeIn plot.bfem Plot the criterion xult 
#' @export
plot_crit = function(x, crit = NULL, ...) {
  color_palette = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", 
    "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
  shape_palette = factor(1:12)
  if(is.null(crit)) crit = x$crit
  df = x$allCriteria
  df$K = as.integer(df$K)
  gg = ggplot2::ggplot(df, 
                       ggplot2::aes_string(x="K", y=crit, shape = "model", 
                                  linetype = "model", col = "model")) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::ylab(base::casefold(crit, upper = T)) +
    ggplot2::scale_shape_manual(values = shape_palette) +
    ggplot2::scale_linetype_manual(values = shape_palette) +
    ggplot2::scale_color_manual(values = color_palette)
  
  gg
}
