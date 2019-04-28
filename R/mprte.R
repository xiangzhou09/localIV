#' Estimation of Marginal Policy Relevant Treatment Effects (MPRTE)
#'
#' \code{mprte} is a function that estimates a class of marginal policy relevant
#' treatment effects (MPRTE) considered in Zhou and Xie (2019). The user needs to
#' specify a \code{policy} as a scalar function of the propensity score.
#'
#' @param object An object of class \code{mte} returned by \code{\link{mte}}.
#' @param policy A univariate scalar function that measures the intensity
#'   of policy intervention across individuals with different levels of
#'   the propensity score.
#'
#' @return Estimate of MPRTE.
#' @export
#'
#' @examples
#' mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   method = "localIV", data = toydata)
#'
#' mprte1 <- mprte(mte_fit, policy = function(p) 1)
#' mprte2 <- mprte(mte_fit, policy = function(p) p)
#' mprte3 <- mprte(mte_fit, policy = function(p) I(p<0.2))
#' c(mprte1, mprte2, mprte3)
#'
#' @references Zhou, Xiang and Yu Xie. 2019. "Marginal Treatment Effects from
#'   A Propensity Score Perspective." Journal of Political Economy.
#'
mprte <- function(object, policy){

  if(!inherits(object, "mte")) stop("object must be an object of class `mte`.")
  if(!is.function(policy)) stop("policy must be a function.")

  # weights
  policy <- Vectorize(policy)
  w <- policy(object$ps)

  # mte_x and mte_u_at_p
  y1_fitted <- as.numeric(object$X[, -1, drop = FALSE] %*% object$coefs$beta1)
  y2_fitted <- as.numeric(object$X[, -1, drop = FALSE] %*% object$coefs$beta2)
  mte_x <- y2_fitted - y1_fitted
  mte_u_at_p <- object$ufun(object$ps)

  # empirical MPRTE
  mprte_emp <- mte_x + mte_u_at_p

  weighted.mean(mprte_emp, w, na.rm = TRUE)
}
