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
  coefs <- list(gamma = gamma, beta1 = beta1, beta2 = beta2,
                theta1 = theta1, theta2 = theta2)

  # MTE mat
  y1_fitted <- as.numeric(X %*% beta1)
  y2_fitted <- as.numeric(X %*% beta2)
  mte_x <- y2_fitted - y1_fitted
  us <- seq(0.005, 0.995, 0.01)
  mte_u <- as.numeric((theta2 - theta1) %*% qnorm(us))
  mte_mat <- outer(mte_x, mte_u, `+`)

  # MTR given x and u
  mtr <- function(x, u, d){
    if(!(d %in% c(0, 1))) stop("d must be 0 or 1")
    if(d == 0) out <- x %*% beta1 + theta1 * qnorm(u) else
      out <- x %*% beta2 + theta2 * qnorm(u)
    as.numeric(out)
  }
  mte <- function(x, u) mtr(x, u, d = 1) - mtr(x, u, d = 0)

  # TR given p
  m1 <- mgcv::gam(y1_fitted ~ s(ps))
  m2 <- mgcv::gam(y2_fitted ~ s(ps))
  pnew <- seq(min(ps), max(ps), length.out = 500)
  pred1 <- mgcv::predict.gam(m1, newdata = list(ps = pnew))
  pred2 <- mgcv::predict.gam(m2, newdata = list(ps = pnew))
  pfun1 <- approxfun(pnew, pred1)
  pfun2 <- approxfun(pnew, pred2)

  # MTR and MTE given p and u
  mtr_tilde <- function(p, u, d){
    if(!(d %in% c(0, 1))) stop("d must be 0 or 1")
    if(d == 0) out <- pfun1(p) + theta1 * qnorm(u) else
      out <- pfun2(p) + theta2 * qnorm(u)
    as.numeric(out)
  }
  mte_tilde <- function(p, u) mtr_tilde(p, u, d = 1) - mtr_tilde(p, u, d = 0)

  # output
  out <- list(mte = mte, mte_tilde = mte_tilde,
              mtr = mtr, mtr_tilde = mtr_tilde,
              coefs = coefs, mte_mat = mte_mat,
              ps = ps, ps_model = normal_switch)
}

