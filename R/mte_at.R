#' Evaluate Marginal Treatment Effects from a Fitted MTE Model.
#'
#' \code{mte_at} evaluates marginal treatment effects at different
#' values of the covariates \code{x} and the latent resistance \code{u}.
#'
#' @param x Values of the pretreatment covariates at which \eqn{\textup{MTE}(x, u)} is evaluated. It
#'   should be a numeric matrix where the number of columns is one less than the number
#'   of columns of the design matrix \eqn{X} in the outcome model. Default is a one-row
#'   matrix where each covariate is set at its sample mean.
#' @param u A numeric vector. Values of the latent resistance \eqn{u} at which
#'   \eqn{\textup{MTE}(x, u)} is evaluated. Note that the estimation involves extrapolation
#'   when the specified u values lie outside of the support of the propensity score.
#' @param model A fitted MTE model returned by \code{\link{mte}}.
#'
#' @return \code{mte_at} returns a data frame. The first \code{length(x)} columns
#'   reflect input values of \code{x}.The following columns are
#'   \item{u}{input values of \code{u}.}
#'   \item{x_comp}{the x-component of the estimated \eqn{\textup{MTE}(x, u)}}
#'   \item{u_comp}{the u-component of the estimated \eqn{\textup{MTE}(x, u)}}
#'   \item{value}{estimated values of \eqn{\textup{MTE}(x, u)}}
#' @export
#' @examples
#' mod <- mte(selection = d ~ x + z, outcome = y ~ x, data = toydata)
#'
#' mte_vals <- mte_at(u = seq(0.005, 0.995, 0.01), model = mod)
#' if(require("ggplot2")){
#'   ggplot(mte_vals, aes(x = u, y = value)) +
#'   geom_line(size = 1) +
#'   xlab("Latent Resistance U") +
#'   ylab("Estimates of MTE at Average values of X") +
#'   theme_minimal(base_size = 14)
#' }
#'
mte_at <- function(x = NULL, u, model){

  # check u and model
  if(missing(u)) stop("`u` must be provided.")
  if(missing(model) || !inherits(model, "mte")){
    stop("model must be an object of class `mte`.")
  }

  # covariate values
  X <- model.matrix(formula(model$mf_o), model$mf_o)
  x <- x %||% matrix(colMeans(X)[-1], nrow = 1)
  if(ncol(x) != ncol(X) - 1){
    stop("`x` must be of length `ncol(X)-1`")
  }
  colnames(x) <- colnames(X)[-1]

  # u component
  mte_u <- model$ufun(u)
  ave_mte_u <- mean(mte_u)
  u_comp <- mte_u - ave_mte_u

  # x component
  mte_x <- as.double(x %*% (model$coefs$beta2 - model$coefs$beta1))
  x_comp <- mte_x + ave_mte_u

  # estimates of MTE given each combination of x and u
  nx <- nrow(x)
  nu <- length(u)
  x_expanded <- x[rep(1:nx, each = length(u)), , drop = FALSE]
  u_expanded <- rep(u, nx)
  x_comp <- rep(x_comp, each = nu)
  u_comp <- rep(u_comp, nx)
  value <- x_comp + u_comp
  out <- data.frame(x_expanded, u = u_expanded, x_comp, u_comp, value)
  out
}
