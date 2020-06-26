#' Evaluate Marginal Treatment Effects Projected onto the Propensity Score
#'
#' \code{mte_tilde_at} evaluates marginal treatment effects
#' projected onto the estimated propensity score. The projection
#' is done via the function \code{\link[mgcv]{gam}}.
#'
#' @param p A numeric vector. Values of the propensity score at which \eqn{\widetilde{\textup{MTE}}(p, u)}
#'   is evaluated.
#' @param u A numeric vector. Values of the latent resistance at which \eqn{\widetilde{\textup{MTE}}(p, u)}
#'   is evaluated.
#' @param model A fitted MTE model returned by \code{\link{mte}}.
#' @param ... Additional parameters passed to \code{\link[mgcv]{gam}}.
#'
#' @return \code{mte_tilde_at} returns a list of two elements:
#'   \item{df}{A data frame containing five columns:\itemize{
#'     \item{\code{p}}{ input values of \code{p}.}
#'     \item{\code{u}}{ input values of \code{u}.}
#'     \item{\code{p_comp}}{ the p-component of the estimated \eqn{\widetilde{\textup{MTE}}(p, u)}}
#'     \item{\code{u_comp}}{ the u-component of the estimated \eqn{\widetilde{\textup{MTE}}(p, u)}}
#'     \item{\code{value}}{ estimated values of \eqn{\widetilde{\textup{MTE}}(p, u)}}
#'     }
#'   }
#'   \item{proj}{Fitted \code{\link[mgcv]{gam}} model for \eqn{E[\mu_1(X)-\mu_0(X)|P(Z)=p]}}
#'
#' @export
#'
#' @examples
#' mod <- mte(selection = d ~ x + z, outcome = y ~ x, data = toydata)
#'
#' u <- p <- seq(0.05, 0.95, 0.1)
#' mte_tilde <- mte_tilde_at(p, u, model = mod)
#'
#' # heatmap showing MTE_tilde(p, u)
#' if(require("ggplot2")){
#' ggplot(mte_tilde$df, aes(x = u, y = p, fill = value)) +
#'   geom_tile() +
#'   scale_fill_gradient(name = expression(widetilde(MTE)(p, u)), low = "yellow", high = "blue") +
#'   xlab("Latent Resistance U") +
#'   ylab("Propensity Score p(Z)") +
#'   theme_minimal(base_size = 14)
#' }
#'
#' mprte_tilde_df <- subset(mte_tilde$df, p == u)
#'
#' # heatmap showing MPRTE_tilde(p)
#' if(require("ggplot2")){
#' ggplot(mprte_tilde_df, aes(x = u, y = p, fill = value)) +
#'   geom_tile() +
#'   scale_fill_gradient(name = expression(widetilde(MPRTE)(p)), low = "yellow", high = "blue") +
#'   xlab("Latent Resistance U") +
#'   ylab("Propensity Score p(Z)") +
#'   theme_minimal(base_size = 14)
#' }
#'
#' # MPRTE_tilde(p) decomposed into the p-component and the u-component
#' if(require(tidyr) && require(dplyr) && require(ggplot2)){
#' mprte_tilde_df %>%
#'   pivot_longer(cols = c(u_comp, p_comp, value)) %>%
#'   mutate(name = recode_factor(name,
#'          `value` = "MPRTE(p)",
#'          `p_comp` = "p(Z) component",
#'          `u_comp` = "U component")) %>%
#'   ggplot(aes(x = p, y = value)) +
#'   geom_line(aes(linetype = name), size = 1) +
#'   scale_linetype(name = "") +
#'   xlab("Propensity Score p(Z)") +
#'   ylab("Treatment Effect") +
#'   theme_minimal(base_size = 14) +
#'   theme(legend.position = "bottom")
#' }
#'
#' @references Zhou, Xiang and Yu Xie. 2019. "\href{https://www.journals.uchicago.edu/doi/abs/10.1086/702172}{Marginal Treatment Effects from A Propensity Score Perspective.}"
#'   Journal of Political Economy, 127(6): 3070-3084.
#' @references Zhou, Xiang and Yu Xie. 2020. "\href{https://journals.sagepub.com/doi/abs/10.1177/0081175019862593}{Heterogeneous Treatment Effects in the Presence of Self-selection:
#'   a Propensity Score Perspective.}" Sociological Methodology.
#'
mte_tilde_at <- function(p, u, model, ...){

  if(missing(model) || !inherits(model, "mte")){
    stop("model must be an object of class `mte`.")
  }

  # u component
  mte_u <- model$ufun(u)
  ave_mte_u <- mean(mte_u)
  u_comp <- mte_u - ave_mte_u

  # p component
  X <- model.matrix(formula(model$mf_o), model$mf_o)
  mte_X <- as.double(X[, -1, drop = FALSE] %*%
                       (model$coefs$beta2 - model$coefs$beta1))
  ps <- model$ps
  proj <- mgcv::gam(mte_X ~ s(ps), ...)
  mte_p <- as.double(mgcv::predict.gam(proj, newdata = list(ps = p)))
  p_comp <- mte_p + ave_mte_u

  # estimates of MTE_tilde given each combination of p and u
  df <- expand.grid(p = p, u = u)
  df$p_comp <- rep(p_comp, length(u))
  df$u_comp <- rep(u_comp, each = length(p))
  df$value <- df$p_comp + df$u_comp

  out <- NULL
  out$df <- df
  out$proj <- proj
  out
}
