# To draw the stream lines as in the textbook, interchange x and z axes. 
# The diagram is drawn in XZ plane.

z <- seq(from = -2, to = 2, len = 1000)
x <- seq(from = 1e-6, to = 2, len = 1000)

f <- function(x, z) {
	x^2/(x^2 + z^2)^(3/2)
}

psi <- outer(x, z, FUN = 'f')

lv <- c(1e-6, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10) # Values of level curves

# To draw the contours with z axis set horizontally, we transpose the matrix
# psi.
contour(t(psi), levels = lv, xlim = c(0.4, 0.6), ylim = c(0, 0.25))
title('Streamlines of a source doublet')

