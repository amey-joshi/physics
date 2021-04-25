x <- seq(from = 0, to = 3, len = 100)
y <- pnorm(x) - pnorm(0)
X <- c(-rev(x), x)
# Recall that we are in the frame of reference such that upper 
# fluid moves with velocity U and the lower one with velocity -U.
Y <- c(-rev(y), y)
plot(
  Y,
  X,
  type = 'l',
  main = 'Transition layer',
  xlab = 'u/U',
  ylab = expression(y / sqrt(4 * nu * t))
)
abline(h = 0, lty = 2)
abline(v = 0, lty = 2)
