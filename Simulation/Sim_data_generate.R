# Data Generation ---------------------------------------------------------

## piece wise responses within a cluster {⇣(s) : s 2 Sj} are identically distributed and have different means across clusters.

#' @  k_true number of cluster 
#' @ para n number of locations 
#' @ nt number of time points
rm(list = ls())

packages = c("wavethresh", "gridExtra","tidyverse","RColorBrewer",'INLA', "deldir", "fields", "igraph","ggraph", "RANN" , "grid", "parallel", "class", "sf", "ggplot2", 'fda')
## Now load or install&load all packages
package.check <- lapply( packages,
                         FUN = function(x) {
                           if (!require(x, character.only = TRUE)) { install.packages(x, dependencies = TRUE) 
                             library(x, character.only = TRUE)
                           }else
                             library(x, character.only = TRUE)
                         } )
library(Matrix.utils)
library(remotes)
remotes::install_github("wmgeolab/rgeoboundaries")
library(rgeoboundaries)
boundary <- geoboundaries("Brazil","adm2")
ggplot(data = boundary) +
  geom_sf()

coords <- st_coordinates(st_centroid(boundary)[,3])
coords = matrix(coords, ncol = 2)
colnames(coords) = c('lon', 'lat')

library(readxl)
brazil_2021 <- read_excel("data/brazil_2021.xls")
# partition and true label---------------------------------------------------------------
ns = nrow(coords)
nt = 256
cluster_true = rep(0, ns)
bbox <- st_bbox(boundary)
center <- c( -45,-19)
idx1 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] >= 0) & (coords[, 2] - center[2] >= 0)
idx2 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] >= 0) & (coords[, 2] - center[2] < 0) 
idx3 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] < 0) & (coords[, 2] - center[2] < 0) 
idx4 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] < 0) & (coords[, 2] - center[2] >= 0) 
idx5 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 < 3.5^2) 
cluster_true[idx1] = 1
a = length(cluster_true[idx1])
cluster_true[idx2] = 2
b = length(cluster_true[idx2])
cluster_true[idx3] = 3
c = length(cluster_true[idx3])
cluster_true[idx4] = 4
d =length(cluster_true[idx4])
cluster_true[idx5]= 5   
e= length(cluster_true[idx5])
c(a,b,c,d,e)
sum(a,b,c,d,e)
k_true = max(cluster_true)

t = seq(0,1, length.out = nt)

J = 8
splinebasis = NULL
sigma2_y <- rep(0.01, k_true)
# splinebasis = create.polygonal.basis(seq(0,1,0.1))
 splinebasis = create.bspline.basis(c(0,1), J)
tPHI= eval.basis(t, splinebasis) 
diag(cov(tPHI))


scal_para <- seq(0.5,2, length.out = k_true)
shape = seq(1,3, length.out = k_true)
lambda <- matrix(0, k_true, J+1)
beta_random <- matrix(1, k_true, J+1)
eta_true <- matrix(1, k_true, J+1)

set.seed(12345)
for (i in 1:k_true){
  lambda[i, ] <- rgamma(J+1, scale = scal_para[i], shape = shape[i])
  print(i)
  V_vec = (1/(lambda[i,]))*sigma2_y[i]
  beta_random[i,] = sapply(V_vec, function(sd,mean) rnorm(1, 0, sd), mean=0 )
}
beta_true <- matrix(1, k_true, J+1)
beta_true[1, ] = c(0, 0.05 * 1:8 - 1) + beta_random[1, ]
beta_true[2, ] = c(0, -0.05 * 1:8 -0.5) + beta_random[2, ]
beta_true[3, ] = c(0, -0.05 * (-4:3)^2) -0.3  + beta_random[3, ]
beta_true[4, ] = c(0, 0.03 * (-4:3)^2 - 1.1) + beta_random[4, ]
 beta_true[5, ] =  beta_random[5, ]-0.4

 mean_true = beta_true %*% t(cbind(1,tPHI))
 
 dim(mean_true)

 eta_true = matrix(0,nrow = ns, ncol = nt)

 eta_true[idx1, ] = rep(mean_true[1, ], nrow(eta_true[idx1, ])) 
 eta_true[idx2, ] = rep(mean_true[2, ], nrow(eta_true[idx2, ]))
 eta_true[idx3, ] = rep(mean_true[3, ], nrow(eta_true[idx3, ]))
 eta_true[idx4, ] = rep(mean_true[4, ], nrow(eta_true[idx4, ]))
 eta_true[idx5, ] = rep(mean_true[5, ], nrow(eta_true[idx5, ]))
 
# # mean_true = beta_true %*% t(tPHI)
# range = range(exp(mean_true))
# plot(t,exp(mean_true[1,]), type = 'l', ylim = range)
# lines(t, exp(mean_true[2,]))
# lines(t,exp(mean_true[3,]))
# lines(t,exp(mean_true[4,]))
# lines(t,exp(mean_true[5,]))
# 
# 
Sigma2 <- matrix(0,nrow = ns, ncol = nt)
Sigma2[idx1, ] <-rnorm(length(Sigma2[idx1]), 0, sigma2_y[1])
Sigma2[idx2, ] <-rnorm(length(Sigma2[idx2]), 0, sigma2_y[2])
Sigma2[idx3, ] <-rnorm(length(Sigma2[idx3]), 0, sigma2_y[3])
Sigma2[idx4, ] <-rnorm(length(Sigma2[idx4]), 0, sigma2_y[4])
Sigma2[idx5, ] <-rnorm(length(Sigma2[idx5]), 0, sigma2_y[5])

population = as.numeric(brazil_2021$`POPULAÇÃO ESTIMADA`)
population = replace_na(population, median(population,na.rm = T))

mu = matrix(0,nrow = ns, ncol = nt)
for(i in 1:ns){mu[i,] = as.matrix(exp(eta_true + Sigma2))[i,]*population[i]}
range(eta_true + Sigma2)
Y = rpois(length(mu),lambda = mu)

Y_obs = matrix(Y, nrow = ns, ncol = nt)

raw_data = cbind(coords, cluster_true,population,Y_obs)

save(raw_data,file = 'data/raw_data_eqvar.Rdata')

ggplot(data = boundary) +
  geom_sf()
