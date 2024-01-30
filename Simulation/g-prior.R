library(parallel)
library(Matrix)
library(statmod)
library(runjags)
library(rjags)
library(INLA)
library(sn)
library(mvtnorm)

######Compute the log-likelihood value for a cluster
######k: the cluster ID;
######Y: observation vector;
######X: inverse discrete wavelet transform matrix;
######lambda: the (J+1) by one vector of scale parameters;
######sigmasq_y: the measurement-error variance;
######ind_gam: the m by one vector of indicator variables;
######membership: the n b one vector of cluster memberships

evaLogLike_Data(Yk,X) <- function
