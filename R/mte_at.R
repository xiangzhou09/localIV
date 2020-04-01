#' Evaluate Marginal Treatment Effects from a Fitted MTE Model.
#'
#' \code{mte_at} evaluates marginal treatment effects at different
#' values of the latent resistance \code{u} with a given \eqn{X=x}.
#'
#' @param x Values of the pretreatment covariates at which \eqn{\textup{MTE}(x, u)} is evaluated. It
#'   should be a numeric vector whose length is one less than the number
#'   of columns of the design matrix \eqn{X} in the outcome model. Default is the sample means.
#' @param u A numeric vector. Values of the latent resistance \eqn{u} at which
#'   \eqn{\textup{MTE}(x, u)} is evaluated. Note that the estimation involves extrapolation
#'   when the specified u values lie outside of the support of the propensity score.
#' @param model A fitted MTE model returned by \code{\link{mte}}.
#'
#' @return \code{mte_at} returns a data frame.
#'   \item{u}{input values of \code{u}.}
#'   \item{x_comp}{the x-component of the estimated \eqn{\textup{MTE}(x, u)}}
#'   \item{u_comp}{the u-component of the estimated \eqn{\textup{MTE}(x, u)}}
#'   \item{value}{estimated values of \eqn{\textup{MTE}(x, u)}}
#' @export
#' @examples
#' mod <- mte(selection = d ~ x + z, outcome = y ~ x, data = toydata)
#'
#' mte_vals <- mte_at(u = seq(0.05, 0.95, 0.1), model = mod)
#' if(require("ggplot2")){
#'   ggplot(mte_vals, aes(x = u, y = value)) +
#'   geom_line(size = 1) +
#'   xlab("Latent Resistance U") +
#'   ylab("Estimates of MTE at Mean Values of X") +
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
  x <- x %||% colMeans(X)[-1]
  if(length(x) != ncol(X) - 1){
    stop("`x` must be of length `ncol(X)-1`")
  }

  # u component
  mte_u <- model$ufun(u)
  ave_mte_u <- mean(mte_u)
  u_comp <- mte_u - ave_mte_u

  # x component
  mte_x <- as.double(x %*% (model$coefs$beta2 - model$coefs$beta1))
  x_comp <- mte_x + ave_mte_u

  # estimates of MTE at different values of u
  value <- x_comp + u_comp
  out <- data.frame(u, x_comp, u_comp, value)
  out
}
