library("splines")
library(INLA)
knots = seq(0,1, length.out  = J)

# estimate mlik for kth cluster using INLA 
evalLogLike_each_INLA(11, Y, tPHI, lambda, sigmasq_y, population, membership)
