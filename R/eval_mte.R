#' Evaluate Marginal Treatment Effects from a Fitted MTE Model.
#'
#' \code{eval_mte} is a function that evaluates marginal treatment effects at any
#' combination of covariates \code{x} and latent resistance \code{u} from a fitted
#' \code{mte} object. Note that the estimation may involve substantial extrapolation
#' when the propensity score has a limited support.
#'
#' @param object An object of class \code{mte} returned by \code{\link{mte}}.
#' @param x A set of pretreatment covariates at which MTE(x, u) is evaluated. Default
#'   is the sample means.
#' @param u Value(s) of the latent resistance \eqn{u} at which MTE(x, u) is evaluated.
#'
#' @return A list of three elements.
#'   \item{mte}{Estimates of MTE(x, u)}
#'   \item{x_comp}{Estimates of \eqn{\mu_1(x)-\mu_0(x)}}
#'   \item{u_comp}{Estimates of \eqn{E[\eta|U=u]}}
#' @export
#'
#' @examples
#' mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   method = "localIV", data = toydata)
#'
#' # plot MTE(x, u) as a function of u
#' u <- seq(0.005, 0.995, 0.01)
#' out <- eval_mte(mte_fit, u = u)
#' plot(out$mte ~ u, type = "l", lwd = 2)
#'
eval_mte <- function(object, x = colMeans(object$X)[-1], u){

  if(!inherits(object, "mte")) stop("object must be an object of class `mte`.")
  if(length(x) != ncol(object$X) - 1) stop("`x` must be of length `ncol(object$X)-1`")

  mte_x <- as.numeric(x %*% (object$coefs$beta2 - object$coefs$beta1))
  mte_u <- object$ufun(u)

  out <- NULL
  out$mte <- mte_x + mte_u
  out$x_comp <- mte_x + mean(mte_u)
  out$u_comp <- mte_u - mean(mte_u)

  out
}
