#' Estimation of Marginal Treatment Effects (MTE)
#'
#' \code{mte} is a function that estimates MTE using either local instrumental
#' variables (local IV) or a normal selection model (Heckman, Urzua, Vytlacil 2006).
#' The user supplies a formula for the treatment selection model, a formula for the
#' outcome model, and a data frame containing the variables. The function returns an
#' object of class \code{mte}. Observations which contain NA (either in \code{selection} or
#' \code{outcome}) are removed.
#'
#' @param selection A formula representing the selection equation.
#' @param outcome A formula representing the outcome equation where the left hand side
#'   is the observed outcome and the right hand side includes predictors of both potential
#'   outcomes.
#' @param data An optional data frame, list, or environment containing the variables
#'   in the model.
#' @param method How to estimate the model: either "\code{localIV}" for local instrumental
#'   variables or "\code{normal}" for a normal selection model.
#' @param bw Bandwidth used for the local polynomial regression in the local IV approach.
#'   Default is 0.25.
#'
#' @return An object of class \code{mte}.
#'  \item{mte}{Fitted MTE function, taking two arguments, \code{x}
#'    for a vector of observed covariates and \code{u} for the latent resistance.}
#'  \item{mte_tilde}{Fitted MTE function, taking two arguments, \code{p} for the propensity
#'    score and \code{u} for the latent resistance. See Zhou and Xie (Forthcoming) for a definition.}
#'  \item{mtr}{Fitted marginal treatment response (MTR) function, taking three arguments, \code{x}
#'    for a vector of observed covariates, \code{u} for the latent resistance, and \code{d} for
#'    treatment status (1 or 2). Available only when \code{method = "normal"}.
#'    See Mogstad, Santos, and Torgovitsky (Forthcoming) for a definition.}
#'  \item{mtr_tilde}{Fitted marginal treatment response (MTR) function, taking three arguments,
#'    \code{p} for the propensity score, \code{u} for the latent resistance, and \code{d} for
#'    treatment status (1 or 2). Available only when \code{method = "normal"}.}
#'  \item{coefs}{A list of fitted coefficients: \code{gamma} for the treatment selection model
#'     (a probit model), \code{beta1} for the baseline outcome, \code{beta2} for the treated outcome,
#'     and \code{theta1} and \code{theta2} for the error covariances when \code{method = "normal"}.}
#'  \item{mte_mat}{An N-by-100 matrix representing estimated MTE for all units, where N
#'    is the sample size. Each row represents a unit and each column represents a value of the
#'    normalized latent resistance, which runs from 0.005 to 0.995 with a step size of 0.01.}
#'  \item{ps}{Estimated propensity scores.}
#'  \item{ps_model}{The propensity score model, an object of class \code{\link[stats]{glm}}
#'     if \code{method = "localIV"}, or an object of class \code{\link[sampleSelection]{selection}}
#'     if \code{method = "normal"}.}
#'  \item{Z}{The model matrix for the treatment selection equation.}
#'  \item{D}{The response vector for the treatment selection equation.}
#'  \item{X}{The model matrix for the outcome equation.}
#'  \item{Y}{The observed outcome.}
#'  \item{call}{The matched call.}
#' @import stats
#' @export
#'
#' @examples
#' mte_fit <- mte(selection = d ~ x + z, outcome = y ~ x,
#'   method = "normal", data = toydata)
#' summary(mte_fit$ps_model)
#'
#' op <- par(mfrow = c(1, 3))
#'
#' # heterogeneous treatment effects by the propensity score
#' with(mte_fit, curve(mte_tilde(p = x, u = 0.5), 0, 1))
#'
#' # heterogeneous treatment effects by the latent resistance
#' with(mte_fit, curve(mte_tilde(p = 0.5, u = x), 0, 1))
#'
#' # heterogeneous treatment effects among marginal entrants
#' with(mte_fit, curve(mte_tilde(p = x, u = x), 0, 1))
#'
#' par(op)
#'
#' @references Heckman, James J., Sergio Urzua, and Edward Vytlacil. 2006.
#'   "Understanding Instrumental Variables in Models with Essential Heterogeneity."
#'   The Review of Economics and Statistics 88:389-432.
#' @references Zhou, Xiang and Yu Xie. Forthcoming. "Marginal Treatment Effects from
#'   A Propensity Score Perspective." Journal of Political Economy.
#' @references Mogstad, Santos, and Torgovitsky. Forthcoming. "Using Instrumental Variables for
#'   Inference About Policy Relevant Treatment Effects." Econometrica.
#'
mte <- function(selection, outcome, data, method = c("localIV", "normal"), bw = 0.25){

  # set data to the parent environment if missing
  if(missing(data)) data <- parent.frame()

  # matched call
  cl <- match.call()

  # model frame for treatment selection model
  m <- match(c("selection", "data"), names(cl), 0)
  mfS <- cl[c(1L, m)]
  mfS[[1L]] <- quote(model.frame)
  names(mfS)[2L] <- "formula"
  mfS$drop.unused.levels <- TRUE
  mfS$na.action <- na.pass
  mfS <- eval(mfS, parent.frame())
  mtS <- terms(mfS)
  Z <- model.matrix(mtS, mfS)
  D <- model.response(mfS)
  DLevels <- levels(as.factor(D))
  D <- as.integer(D == utils::tail(DLevels, 1L))

  # model frame for outcome model
  m <- match(c("outcome", "data"), names(cl), 0)
  mfO <- cl[c(1L, m)]
  mfO[[1L]] <- quote(model.frame)
  names(mfO)[2L] <- "formula"
  mfO$drop.unused.levels <- TRUE
  mfO$na.action <- na.pass
  mfO <- eval(mfO, parent.frame())
  mtO <- terms(mfO)
  X <- model.matrix(mtO, mfO)
  Y <- model.response(mfO)

  # delete rows with missing data
  badRow <- is.na(D) | is.infinite(D) | is.na(Y) | is.infinite(Y)
  badRow <- badRow | apply(Z, 1, function(v) any(is.na(v) | is.infinite(v)))
  badRow <- badRow | apply(X, 1, function(v) any(is.na(v) | is.infinite(v)))

  Z <- Z[!badRow, , drop = FALSE]
  D <- D[!badRow]
  X <- X[!badRow, , drop = FALSE]
  Y <- Y[!badRow]
  N <- length(D)

  environment(mte_normal) <- environment()

  method <- match.arg(method, c("localIV", "normal"))
  if (method == "normal"){
    out <- mte_normal(selection, outcome, data)
  } else {
    out <- mte_localIV(Z, D, X, Y, bw = bw)
  }
  out$Z <- Z
  out$D <- D
  out$X <- X
  out$Y <- Y
  out$cl <- cl
  class(out) <- c("mte", "list")
  out
}




