#' Evaluate Marginal Treatment Effects Projected onto the Propensity Score
#'
#' \code{eval_mte_tilde} is a function that evaluates marginal treatment effects
#' projected onto the estimated propensity score (Zhou and Xie 2019). The projection
#' is done via the function \code{\link[mgcv]{gam}} with default parameters.
#'
#' @param object An object of class \code{mte} returned by \code{\link{mte}}.
#' @param p Value(s) of the propensity score \eqn{p} at which MTE_tilde(p, u) is evaluated.
#' @param u Value(s) of the latent resistance \eqn{u} at which MTE_tilde(p, u) is evaluated.
#'
#' @return A list of four elements.
#'   \item{mte_tilde}{Estimates of MTE_tilde(p, u)}
#'   \item{p_comp}{Estimates of \eqn{E[\mu_1(X)-\mu_0(X)|P(Z)=p]}}
#'   \item{u_comp}{Estimates of \eqn{E[\eta|U=u]}}
#'   \item{model}{Fitted model for \eqn{E[\mu_1(X)-\mu_0(X)|P(Z)=p]}}
#' @export
#'
#' @examples
#'  mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'  method = "localIV", data = toydata)
#'
#'  # heatmap showing MTE_tilde(p, u)
#'  library(plotly)
#'  p <- rep(seq(0.05, 0.95, 0.1), 10)
#'  u <- rep(seq(0.05, 0.95, 0.1), each = 10)
#'  out1 <- eval_mte_tilde(mte_fit, p = p, u = u)
#'  plot_ly(x = u, y = p, z = out1$mte_tilde, type = "heatmap")
#'
#'  # heatmap showing MPRTE_tilde(p)
#'  p <- seq(0.05, 0.95, 0.1)
#'  u <- p
#'  out2 <- eval_mte_tilde(mte_fit, p = p, u = u)
#'  plot_ly(x = u, y = p, z = out2$mte_tilde, type = "heatmap")
#'
#'  # decompose MPRTE_tilde(p) into the p-component and the u-component
#'  y <- with(out2, cbind(mte_tilde, p_comp, u_comp))
#'  matplot(x = p, y = y, type = "l", lwd = 2)
#'
#' @references Zhou, Xiang and Yu Xie. 2019. "Marginal Treatment Effects from
#'   A Propensity Score Perspective." Journal of Political Economy.
#'
eval_mte_tilde <- function(object, p, u){

  if(!inherits(object, "mte")) stop("object must be an object of class `mte`.")

  X <- object$X[, -1, drop = FALSE]
  mte_X <- as.numeric(X %*% (object$coefs$beta2 - object$coefs$beta1))
  ps <- object$ps
  proj <- mgcv::gam(mte_X ~ s(ps))

  mte_p <- as.numeric(mgcv::predict.gam(proj, newdata = list(ps = p)))
  mte_u <- object$ufun(u)

  out <- NULL
  out$mte_tilde <- mte_p + mte_u
  out$p_comp <- mte_p + mean(mte_u)
  out$u_comp <- mte_u - mean(mte_u)
  out$model <- proj

  out
}
