#' @rdname mte
#' @export
mte_normal <- function(mf_s, mf_o){

  # extract N, Z, D, X and Y
  N <- nrow(mf_s)
  Z <- model.matrix(formula(mf_s), mf_s)
  D <- model.response(mf_s)
  if(length(unique(D))!=2) stop("`D` must be binary.")
  X <- model.matrix(formula(mf_o), mf_o)
  Y <- model.response(mf_o)
  pZ <- ncol(Z)
  pX <- ncol(X)

  # formulas
  selection <- formula(mf_s)
  outcome <- formula(mf_o)

  # data
  names_in_s_only <- setdiff(names(mf_s), names(mf_o))
  data <- cbind(mf_s[names_in_s_only], mf_o)

  # initialization via Heckman two-step estimation
  twoStep <- sampleSelection::heckit5fit(selection, outcome, outcome, data=data)
  ind <- twoStep$param$index
  start <- coef(twoStep, part="full")[c(ind$betaS, ind$betaO1, ind$sigma1,
                                        ind$rho1, ind$betaO2, ind$sigma2, ind$rho2)]
  names(start) <- sub( "^[SO][12]?:", "", names(start))

  # switching regression model
  YO1 <- YO2 <- Y
  XO1 <- XO2 <- X
  YO1[D == 1] <- NA
  YO2[D == 0] <- NA
  XO1[D == 1, ] <- NA
  XO2[D == 0, ] <- NA
  normal_switch <- sampleSelection::tobit5fit(D, Z, YO1, XO1, YO2, XO2, start=start)

  # propensity score from the normal selection model
  gamma <- normal_switch$estimate[1:pZ]
  ps <- as.double(pnorm(Z %*% matrix(gamma)))

  # parametric MTE estimation
  beta1 <- normal_switch$estimate[(pZ + 1):(pZ + pX)]
  beta2 <- normal_switch$estimate[(pZ + pX + 3):(pZ + 2 * pX + 2)]
  sigma1 <- normal_switch$estimate["sigma1"]
  rho1 <- normal_switch$estimate["rho1"]
  theta1 <- c(- sigma1 * rho1, use.names = FALSE)
  sigma2 <- normal_switch$estimate["sigma2"]
  rho2 <- normal_switch$estimate["rho2"]
  theta2 <- c(- sigma2 * rho2, use.names = FALSE)

  # all coefficients
  coefs <- list(gamma = gamma,
                beta10 = beta1[1],
                beta1 = beta1[-1],
                beta20 = beta2[1],
                beta2 = beta2[-1],
                theta1 = theta1,
                theta2 = theta2)

  # u component of mte(x, u)
  ufun <- function(u) beta20 - beta10 + (theta2 - theta1) * qnorm(u)
  fn_env(ufun) <- new_environment(coefs, parent = global_env())

  # output
  out <- list(coefs = coefs, ufun = ufun, ps = ps, ps_model = normal_switch,
              mf_s = mf_s, mf_o = mf_o, selection = selection, outcome = outcome)
  class(out) <- "mte"
  out
}
