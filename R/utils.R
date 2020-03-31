utils::globalVariables(c("beta10", "beta20", "theta1", "theta2"))

# logical or infix function
`%||%` <- function(a, b) if (!is.null(a)) a else b
