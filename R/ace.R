#' Estimating Average Causal Effects from a Fitted MTE Model.
#'
#' \code{ace} estimates Average Causal Effects (ACE) from a fitted MTE model.
#' The estimand can be average treatment effect (ATE), average treatment effect on the treated (ATT),
#' average treatment effect on the untreated (ATU), or the Marginal Policy Relevant
#' Treatment Effect (MPRTE) defined in Zhou and Xie (2019).
#'
#' @param model A fitted \code{mte} model returned by \code{\link{mte}}.
#' @param estimand Type of estimand: \code{"ate"}, \code{"att"}, \code{"atu"}, or \code{"mprte"}.
#' @param policy An \code{\link{expression}} written as a function of \code{p}. This is used
#'   only when \code{estimand="mprte"}.
#'
#' @return Estimate of ATE, ATT, ATU, or MPRTE
#' @export
#'
#' @examples
#' mod <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   data = toydata)
#'
#' ate <- ace(mod, "ate")
#' att <- ace(mod, "att")
#' atu <- ace(mod, "atu")
#' mprte1 <- ace(mod, "mprte")
#' mprte2 <- ace(mod, "mprte", policy = p)
#' mprte3 <- ace(mod, "mprte", policy = 1-p)
#' mprte4 <- ace(mod, "mprte", policy = I(p<0.25))
#' c(ate, att, atu, mprte1, mprte2, mprte3, mprte4)
#'
#' @references Heckman, James J., Sergio Urzua, and Edward Vytlacil. 2006.
#'   "Understanding Instrumental Variables in Models with Essential Heterogeneity."
#'   The Review of Economics and Statistics 88:389-432.
#' @references Zhou, Xiang and Yu Xie. 2019. "\href{https://www.journals.uchicago.edu/doi/abs/10.1086/702172}{Marginal Treatment Effects from A Propensity Score Perspective.}"
#'   Journal of Political Economy, 127(6): 3070-3084.
#' @references Zhou, Xiang and Yu Xie. 2020. "\href{https://scholar.harvard.edu/files/xzhou/files/zhou-xie2019_hte.pdf}{Heterogeneous Treatment Effects in the Presence of Self-selection:
#'   a Propensity Score Perspective.}" Sociological Methodology.
ace <- function(model,
                estimand = c("ate", "att", "atu", "mprte"),
                policy = 1){

  if(!inherits(model, "mte")) stop("`model` must be an object of class `mte`.")
  estimand <- match.arg(estimand)

  X <- model.matrix(formula(model$mf_o), model$mf_o)

  # us and MTE mat
  y1_fitted <- as.double(X[, -1, drop = FALSE] %*% model$coefs$beta1)
  y2_fitted <- as.double(X[, -1, drop = FALSE] %*% model$coefs$beta2)
  mte_x <- y2_fitted - y1_fitted

  if(estimand == "mprte"){
    policy_expr <- enexpr(policy)
    policy <- new_function(exprs(p = ), policy_expr)
    w <- Vectorize(policy)(model$ps)
    if(any(w<0)) stop("`policy` should not imply negative weights.")
    if(all(w==0)) stop("`policy` has no empirical support.")
    mte_u_at_p <- model$ufun(model$ps)
    mprte_emp <- mte_x + mte_u_at_p
    out <- weighted.mean(mprte_emp, w, na.rm = TRUE)
    names(out) <- paste0("mprte: ", as_label(policy_expr))
    return(out)
  }

  us <- seq(0.005, 0.995, 0.01)
  mte_u <- model$ufun(us)
  mte_mat <- outer(mte_x, mte_u, `+`)
  D <- model.response(model$mf_s)

  if (estimand == "ate"){
    ate_tilde_p <- rowMeans(mte_mat, na.rm = TRUE)
    out <- mean(ate_tilde_p, na.rm = TRUE)
    names(out) <- "ate"
  } else if (estimand == "att"){
    wtt_tilde <- outer(model$ps, us, `>=`)
    att_tilde_p <- rowSums(mte_mat * wtt_tilde, na.rm = TRUE)/
      rowSums(wtt_tilde, na.rm = TRUE)
    out <- mean(att_tilde_p[D==1], na.rm = TRUE)
    names(out) <- "att"
  } else{
    wtu_tilde <- outer(model$ps, us, `<`)
    atu_tilde_p <- rowSums(mte_mat * wtu_tilde, na.rm = TRUE)/
      rowSums(wtu_tilde, na.rm = TRUE)
    out <- mean(atu_tilde_p[D==0], na.rm = TRUE)
    names(out) <- "atu"
  }
  out
}
