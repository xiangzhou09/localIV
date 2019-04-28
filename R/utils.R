utils::globalVariables(c("Z", "D", "X", "Y", "N"))

mte_normal <- function(selection, outcome, data){

  # normal selection model
  normal_switch <- sampleSelection::selection(selection, list(outcome, outcome),
                                              data, type = 5)

  # propensity score from the normal selection model
  gamma <- normal_switch$estimate[1:ncol(Z)]
  ps <- as.numeric(pnorm(Z %*% matrix(gamma)))

  # parametric MTE estimation
  beta1 <- normal_switch$estimate[(ncol(Z) + 1):(ncol(Z) + ncol(X))]
  beta2 <- normal_switch$estimate[(ncol(Z) + ncol(X) + 3):(ncol(Z) + 2 * ncol(X) + 2)]

  sigma1 <- normal_switch$estimate["sigma1"]
  rho1 <- normal_switch$estimate["rho1"]
  theta1 <- - sigma1 * rho1
  sigma2 <- normal_switch$estimate["sigma2"]
  rho2 <- normal_switch$estimate["rho2"]
  theta2 <- - sigma2 * rho2

  # all coefficients
  coefs <- list(gamma = gamma, beta1 = beta1[-1], beta2 = beta2[-1],
                theta1 = theta1, theta2 = theta2)

  ufun <- function(u) beta2[1] - beta1[1] + (theta2 - theta1) * qnorm(u)

  # output
  out <- list(coefs = coefs, ufun = ufun, ps = ps, ps_model = normal_switch)
}

mte_localIV <- function(selection, data, bw = 0.25){

  # transform X and Z into matrices if needed
  if(!is.matrix(X)) X <- as.matrix(X)
  if(!is.matrix(Z)) Z <- as.matrix(Z)

  # X without constant term and us
  X <- X[, -1, drop = FALSE]
  us <- seq(0.005, 0.995, 0.01)

  # propensity score model
  ps_logit <- glm(formula = selection, family = binomial("probit"), data = data)
  ps <- ps_logit$fitted.values

  # local IV estimation of MTE #

  # step 1: fit Y, X, Xp on p using local linear regressions
  Xp <- X * matrix(rep(ps, ncol(X)), N, ncol(X))
  Y_res <- loess(Y ~ ps, degree = 1)$res
  X_res <- apply(X, 2, function(x) loess(x ~ ps, degree = 1)$res)
  Xp_res <- apply(Xp, 2, function(x) loess(x ~ ps, degree = 1)$res)

  # step 2: double residual regression
  res_reg <- lm.fit(x = cbind(X_res, Xp_res), y = Y_res)
  beta1 <- res_reg$coef[1:ncol(X)]
  beta2 <- res_reg$coef[(ncol(X) + 1):(2 * ncol(X))] + beta1

  # step 3: fit delta and K(p)
  Y_unobserved <- Y - as.numeric(cbind(X, Xp) %*% res_reg$coef)
  mte_u_poly <- KernSmooth::locpoly(ps, Y_unobserved, drv = 1L, bandwidth = bw,
                                    gridsize = 100L, range.x = c(0.005, 0.995))
  # mte_u <- mte_u_poly$y
  ufun <- approxfun(mte_u_poly)

  # all coefficients
  coefs <- list(gamma = ps_logit$coefficients,
                beta1 = beta1, beta2 = beta2,
                theta1 = NULL, theta2 = NULL)
  # output
  out <- list(coefs = coefs, ufun = ufun, ps = ps, ps_model = ps_logit)
}

