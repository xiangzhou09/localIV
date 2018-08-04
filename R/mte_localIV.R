mte_localIV <- function(Z, D, X, Y, bw = 0.25){

  # transform X and Z into matrices if needed
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(Z)) Z <- as.matrix(Z)

  # sample size
  N <- length(D)
  X <- X[, -1, drop = FALSE]
  us <- seq(0.005, 0.995, 0.01)

  # propensity score model
  ps_logit <- glm(D ~ 0 + Z, family = binomial("probit"))
  ps <- ps_logit$fitted.values

  # local IV estimation of MTE #

  # step 1: fit Y, X, Xp on p using local linear regressions
  Xp <- X * matrix(rep(ps, ncol(X)), N, ncol(X))
  Y_res <- loess(Y ~ ps, degree=1)$res
  X_res <- apply(X, 2, function(x) loess(x ~ ps, degree=1)$res)
  Xp_res <- apply(Xp, 2, function(x) loess(x ~ ps, degree=1)$res)

  # step 2: double residual regression
  res_reg <- lm.fit(x = cbind(X_res, Xp_res), y = Y_res)
  beta1 <- res_reg$coef[1:ncol(X)]
  beta2 <- res_reg$coef[(ncol(X)+1):(2*ncol(X))] + beta1
  y1_fitted <- as.numeric(X %*% beta1)
  y2_fitted <- as.numeric(X %*% beta2)
  mte_x <- y2_fitted - y1_fitted

  # step 3: fit delta and K(p)
  Y_unobserved <- Y - as.numeric(cbind(X, Xp) %*% res_reg$coef)
  mte_u_poly <- KernSmooth::locpoly(ps, Y_unobserved, drv = 1L, bandwidth = bw,
                                     gridsize = 100L, range.x = c(0.005, 0.995))
  mte_u <- mte_u_poly$y
  ufun <- approxfun(mte_u_poly)

  # all coefficients
  coefs <- list(gamma = ps_logit$coefficients,
                beta1 = beta1, beta2 = beta2,
                theta1 = NULL, theta2 = NULL)

  # MTE mat
  mte_mat <- outer(mte_x, mte_u, `+`)

  # MTE given x and u
  mte <- function(x, u) as.numeric(x %*% (beta2 - beta1) + ufun(u))

  # TR given p
  m1 <- mgcv::gam(y1_fitted ~ s(ps))
  m2 <- mgcv::gam(y2_fitted ~ s(ps))
  pnew <- seq(min(ps), max(ps), length.out = 500)
  pred1 <- mgcv::predict.gam(m1, newdata = list(ps = pnew))
  pred2 <- mgcv::predict.gam(m2, newdata = list(ps = pnew))
  pfun1 <- approxfun(pnew, pred1)
  pfun2 <- approxfun(pnew, pred2)

  # MTE given p and u
  mte_tilde <- function(p, u) pfun2(p) - pfun1(p) + ufun(u)

  out <- list(mte = mte, mte_tilde = mte_tilde,
              mtr = NULL, mtr_tilde = NULL,
              coefs = coefs, mte_mat = mte_mat,
              ps = ps, ps_model = ps_logit)
}

