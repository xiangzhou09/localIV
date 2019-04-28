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
#' @return A list of two elements.
#'   \item{mte_tilde}{Estimates of MTE_tilde(p, u)}
#'   \item{model}{Fitted model of \eqn{x'(\beta_1 - \beta_0)} as a function of the
#'   propensity score}
#' @export
#'
#' @examples
#' mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   method = "localIV", data = toydata)
#'
#' x <- seq(0.05, 0.95, 0.05)
#' mte_tilde_p <- eval_mte_tilde(mte_fit, p = x, u = 0.5)$mte_tilde
#' mte_tilde_u <- eval_mte_tilde(mte_fit, p = 0.5, u = x)$mte_tilde
#' mprte_tilde_p <- eval_mte_tilde(mte_fit, p = x, u = x)$mte_tilde
#'
#' out <- cbind(mte_tilde_p, mte_tilde_u, mprte_tilde_p)
#' matplot(x = x, y = out, type = "l")
#'
#' @references Zhou, Xiang and Yu Xie. 2019. "Marginal Treatment Effects from
#'   A Propensity Score Perspective." Journal of Political Economy.
eval_mte_tilde <- function(object, p, u){

  if(!inherits(object, "mte")) stop("object must be an object of class `mte`.")
  # if(length(p) != length(u)) stop("`p` and `u` must be of the same length")

  X <- object$X[, -1, drop = FALSE]
  mte_X <- as.numeric(X %*% (object$coefs$beta2 - object$coefs$beta1))
  ps <- object$ps

  proj <- mgcv::gam(mte_X ~ s(ps))
  mte_p <- mgcv::predict.gam(proj, newdata = list(ps = p))

  mte_u <- object$ufun(u)

  list(mte_tilde = as.numeric(mte_p) + mte_u, model = proj)
}
