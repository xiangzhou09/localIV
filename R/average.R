#' Estimation of Average Causal Effects from Marginal Treatment Effects
#'
#' \code{average} is a function that estimates conventional causal parameters including
#' average treatment effect (ATE), average treatment effect on the treated (ATT), and
#' average treatment effect on the untreated (ATU). Note that the estimation may involve
#' substantial extrapolation when the propensity score has a limited support.
#'
#' @param object An object of class \code{mte} returned by \code{\link{mte}}.
#' @param estimand Type of estimand: \code{"ate"}, \code{"att"}, or \code{"atu"}.
#'
#' @return Estimate of ATE, ATT, or ATU.
#' @export
#'
#' @examples
#' mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   method = "localIV", data = toydata)
#'
#' ate <- average(mte_fit, "ate")
#' att <- average(mte_fit, "att")
#' atu <- average(mte_fit, "atu")
#' c(ate, att, atu)
#'
#' @references Heckman, James J., Sergio Urzua, and Edward Vytlacil. 2006.
#'   "Understanding Instrumental Variables in Models with Essential Heterogeneity."
#'   The Review of Economics and Statistics 88:389-432.

average <- function(object, estimand = c("ate", "att", "atu")){

  if(!inherits(object, "mte")) stop("object must be an object of class `mte`.")
  estimand <- match.arg(estimand, c("ate", "att", "atu"))

  # us and MTE mat
  y1_fitted <- as.numeric(object$X[, -1, drop = FALSE] %*% object$coefs$beta1)
  y2_fitted <- as.numeric(object$X[, -1, drop = FALSE] %*% object$coefs$beta2)
  mte_x <- y2_fitted - y1_fitted
  us <- seq(0.005, 0.995, 0.01)
  mte_u <- object$ufun(us)
  mte_mat <- outer(mte_x, mte_u, `+`)

  # evaluate ate, att, and atu
  if (estimand == "ate"){
    ate_tilde_p <- rowMeans(mte_mat, na.rm = TRUE)
    out <- mean(ate_tilde_p, na.rm = TRUE)
  } else if (estimand == "att"){
    wtt_tilde <- outer(object$ps, us, `>=`)
    att_tilde_p <- rowSums(mte_mat * wtt_tilde, na.rm = TRUE)/
      rowSums(wtt_tilde, na.rm = TRUE)
    out <- mean(att_tilde_p[object$D==1], na.rm = TRUE)
  } else{
    wtu_tilde <- outer(object$ps, us, `<`)
    atu_tilde_p <- rowSums(mte_mat * wtu_tilde, na.rm = TRUE)/
      rowSums(wtu_tilde, na.rm = TRUE)
    out <- mean(atu_tilde_p[object$D==0], na.rm = TRUE)
  }
  out
}
