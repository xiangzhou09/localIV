
# set.seed(1)
# n <- 10000
# vc <- diag(3)
# vc[lower.tri(vc)] <- c(0.2, 0.6, 0.1)
# vc[upper.tri(vc)] <- vc[lower.tri(vc)]
# eps <- MASS::mvrnorm(n, c(0,0,0), vc)
# x <- rnorm(n)
# z <- rnorm(n)
# d <- (x + z + eps[, 1] > 0)
# y0 <- x + eps[, 2]
# y1 <- 2 + 1.5 * x + eps[, 3]
# y <- y1 * d + y0 * (1-d)
# toydata <- data.frame(y = y, x = x, z = z, d = d)
# use_data(toydata)
