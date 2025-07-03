rm(list=ls())

schaffer_f6 <- function(x, y) {
	r <- x**2 + y**2
	a <- sin(sqrt(r))
	b <- 1.0 + 0.001 * r
	return( 0.5 + ( a**2- 0.5 ) / ( b**2 ) )
}

schwefel <- Vectorize(function(x, y) {
	s <- 0.0
	for(xi in c(x,y)) {
		if(xi > 500 || xi < -500) {
			s <- s + 0.02 * xi * xi
		}
		else {
			s <- s + -xi * sin(sqrt(abs(xi)))
		}
	}
	return( 418.9829 * 2 + s )
})

x <- seq(-500, 500, 5)
y <- seq(-500, 500, 5)

z <- outer(x, y, schwefel)

persp(x, y, z, 
	xlab='x', ylab='y', 
	zlab='schaffer_f6', 
	main='3D Plot', col='pink', 
	shade=.4,
	theta = 30, phi = 0, ticktype='detailed',
	scale=TRUE,
	expand=0.4)
