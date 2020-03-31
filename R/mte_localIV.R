#' @rdname mte
#' @export
mte_localIV <- function(mf_s, mf_o, bw = NULL){

  # setting default bandwidth
  bw <- bw %||% 0.25

  # extract Z, D, X and Y
  N <- nrow(mf_s)
  Z <- model.matrix(formula(mf_s), mf_s)
  D <- model.response(mf_s)
  if(length(unique(D))!=2) stop("`D` must be binary.")
  X <- model.matrix(formula(mf_o), mf_o)
  Y <- model.response(mf_o)

  # X without constant term and us
  X <- X[, -1, drop = FALSE]
  us <- seq(0.005, 0.995, 0.01)

  # formulas
  selection <- formula(mf_s)
  outcome <- formula(mf_o)

  # data
  names_in_s_only <- setdiff(names(mf_s), names(mf_o))
  data <- cbind(mf_s[names_in_s_only], mf_o)

  # propensity score model
  ps_logit <- glm(selection, family = binomial("probit"), data = mf_s)
  ps <- ps_logit$fitted.values

  # local IV estimation of MTE #
  pX <- ncol(X)

  # step 1: fit Y, X, Xp on p using local linear regressions
  Xp <- X * matrix(rep(ps, pX), N, pX)
  Y_res <- loess(Y ~ ps, degree = 1)$res
  X_res <- apply(X, 2, function(x) loess(x ~ ps, degree = 1)$res)
  Xp_res <- apply(Xp, 2, function(x) loess(x ~ ps, degree = 1)$res)

  # step 2: double residual regression
  res_reg <- lm.fit(x = cbind(X_res, Xp_res), y = Y_res)
  beta1 <- res_reg$coef[1:pX]
  beta2 <- res_reg$coef[(pX + 1):(2 * pX)] + beta1

  # step 3: fit delta and K(p)
  Y_unobserved <- Y - as.double(cbind(X, Xp) %*% res_reg$coef)
  mte_u_poly <- KernSmooth::locpoly(ps, Y_unobserved, drv = 1L, bandwidth = bw,
                                    gridsize = 100L, range.x = c(0.005, 0.995))
  # mte_u <- mte_u_poly$y
  ufun <- approxfun(mte_u_poly)

  # all coefficients
  coefs <- list(gamma = ps_logit$coefficients,
                beta10 = NULL,
                beta1 = beta1,
                beta20 = NULL,
                beta2 = beta2,
                theta1 = NULL,
                theta2 = NULL)

  # output
  out <- list(coefs = coefs, ufun = ufun, ps = ps, ps_model = ps_logit,
              mf_s = mf_s, mf_o = mf_o, selection = selection, outcome = outcome)
  class(out) <- "mte"
  out
}
