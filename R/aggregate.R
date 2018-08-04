#' Estimation of Average Causal Effects from Marginal Treatment Effects
#'
#' \code{average} is a function that estimates conventional causal parameters including
#' average treatment effect (ATE), average treatment effect on the treated (ATT), and
#' average treatment effect on the untreated (ATU). Note that the estimation involves
#' substantial extrapolation when the propensity score has a limited support.
#'
#' @param mte_fit An object of class \code{mte} returned by \code{\link{mte}}.
#' @param estimand Type of estimand: \code{"ate"}, \code{"att"}, or \code{"atu"}.
#'
#' @return Estimate of ATE, ATT, or ATU.
#' @export
#'
#' @examples
#' mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   method = "localIV", data = toydata)
#'
#' ate <- average(mte_fit, estimand = "ate")
#' att <- average(mte_fit, estimand = "att")
#' c(ate, att)
#'
#' @references Heckman, James J., Sergio Urzua, and Edward Vytlacil. 2006.
#'   "Understanding Instrumental Variables in Models with Essential Heterogeneity."
#'   The Review of Economics and Statistics 88:389-432.

average <- function(mte_fit, estimand = c("ate", "att", "atu")){

  if(!inherits(mte_fit, "mte")) stop("mte_fit must be an object of class `mte`.")
  estimand <- match.arg(estimand, c("ate", "att", "atu"))
  us <- get("us", environment(mte_fit$mte))

  if (estimand == "ate"){
    ate_tilde_p <- rowMeans(mte_fit$mte_mat, na.rm = TRUE)
    out <- mean(ate_tilde_p, na.rm = TRUE)
  } else if (estimand == "att"){
    wtt_tilde <- outer(mte_fit$ps, us, `>=`)
    att_tilde_p <- rowSums(mte_fit$mte_mat * wtt_tilde, na.rm = TRUE)/
      rowSums(wtt_tilde, na.rm = TRUE)
    out <- mean(att_tilde_p[mte_fit$D==1], na.rm = TRUE)
  } else{
    wtu_tilde <- outer(mte_fit$ps, us, `<`)
    atu_tilde_p <- rowSums(mte_fit$mte_mat * wtu_tilde, na.rm = TRUE)/
      rowSums(wtu_tilde, na.rm = TRUE)
    out <- mean(atu_tilde_p[mte_fit$D==0], na.rm = TRUE)
  }
  out
}

#' Estimation of Marginal Policy Relevant Treatment Effects (MPRTE)
#'
#' \code{mprte} is a function that estimates a class of marginal policy relevant
#' treatment effects (MPRTE) considered in Zhou and Xie (2018). The user needs to
#' specify a \code{policy} as a scalar function of the propensity score.
#'
#' @param mte_fit An object of class \code{mte} returned by \code{\link{mte}}.
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
#' @references Zhou, Xiang and Yu Xie. Forthcoming. "Marginal Treatment Effects from
#'   A Propensity Score Perspective." Journal of Political Economy.
mprte <- function(mte_fit, policy){
  if(!inherits(mte_fit, "mte")) stop("mte_fit must be an object of class `mte`.")
  if(!is.function(policy)) stop("policy must be a function.")
  policy <- Vectorize(policy)
  x <- mte_fit$mte_tilde(mte_fit$ps, mte_fit$ps)
  w <- policy(mte_fit$ps)
  weighted.mean(x, w, na.rm = TRUE)
}
