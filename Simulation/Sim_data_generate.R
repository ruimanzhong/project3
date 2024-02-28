# Data Generation ---------------------------------------------------------

## piece wise responses within a cluster {⇣(s) : s 2 Sj} are identically distributed and have different means across clusters.

#' @  k_true number of cluster 
#' @ para n number of locations 
#' @ nt number of time points
rm(list = ls())

packages = c("wavethresh", "gridExtra","tidyverse","RColorBrewer",'INLA', "deldir", "fields", "igraph","ggraph", "RANN" , "grid", "parallel", "class", "sf", "ggplot2")
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
# remotes::install_github("wmgeolab/rgeoboundaries")
library(rgeoboundaries)

# State_cluster -----------------------------------------------------------

boundary <-  geoboundaries("Brazil","adm1") 
library(readr)
brazil_pop <- read_csv("~/Documents/project3/data/brazil_pop.csv")
boundary$population = brazil_pop$BrazilStatesPopulation2021


# County_cluster ----------------------------------------------------------
# boundary <- geoboundaries("Brazil","adm2")
# library(readxl)
# brazil_2021 <- read_excel("~/Documents/project3/data/brazil_2021.xls")
# population = as.numeric(brazil_2021$`POPULAÇÃO ESTIMADA`)
# population = replace_na(population, median(population,na.rm = T))
# boundary$population <- population
# boundary <- boundary %>% select(population)


# Data --------------------------------------------------------------------
coords <- st_coordinates(st_centroid(boundary)[,3])
coords = matrix(coords, ncol = 2)
colnames(coords) = c('lon', 'lat')
# partition and true label---------------------------------------------------------------
ns = nrow(coords)
nt = 128
cluster_true = rep(0, ns)
bbox <- st_bbox(boundary)
center <- c( -45,-19)
#idx1 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 < 3.5^2) 
# idx5 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] >= 0) & (coords[, 2] - center[2] >= 0)
#  idx2 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] >= 0) & (coords[, 2] - center[2] < 0) 
# idx3 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] < 0) & (coords[, 2] - center[2] < 0) 
# idx4 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 3.5^2) & (coords[, 1] - center[1] < 0) & (coords[, 2] - center[2] >= 0) 
 idx1 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 < 10^2) 
 idx2 = ((coords[, 1] - center[1])^2 + (coords[, 2] - center[2])^2 >= 10^2) 
cluster_true[idx1] = 1
a = length(cluster_true[idx1])
# cluster_true[idx2] = 2
# b = length(cluster_true[idx2])
# cluster_true[idx3] = 3
# c = length(cluster_true[idx3])
# cluster_true[idx4] = 4
# d =length(cluster_true[idx4])
# cluster_true[idx5]= 5   
# e= length(cluster_true[idx5])
# c(a,b,c,d,e)
# sum(a,b,c,d,e)
cluster_true[idx2]= 2

k_true = max(cluster_true)

t = seq(0,1, length.out = nt)

# J = 8
# splinebasis = NULL
sigma2_y <- rep(0.001, k_true)
# # splinebasis = create.polygonal.basis(seq(0,1,0.1))
#  splinebasis = create.bspline.basis(c(0,1), J)
# tPHI= eval.basis(t, splinebasis) 


# polynomial basis --------------------------------------------------------
J = 3
dim(tPHI)
tPHI = matrix(c(rep(1,nt), t, t^2, t^3), nrow = nt)

# # scal_para <- seq(0.5,2, length.out = k_true)
# shape = seq(1,3, length.out = k_true)
lambda <- matrix(1, k_true, J+1)
eta <- matrix(1, k_true, J+1)


set.seed(12345)
beta_1 = mvrnorm(a, mu = c(-4,rep(1,J)), Sigma = sigma2_y[1]*(diag(1/lambda[1,])) )
beta_2 = mvrnorm(ns-a, mu = rep(-1,J+1), Sigma = sigma2_y[2]*(diag(1/lambda[2,])) )

# all_locations have the same beta ----------------------------------------
# eta = matrix(0,nrow = ns, ncol = nt)
# eta[idx1, ] = t(matrix(rep(beta_1 %*% t(tPHI), nrow(eta_true[idx1, ])), ncol = a))
# eta[idx2, ] = t(matrix(rep(beta_2 %*% t(tPHI), nrow(eta_true[idx2, ])), ncol = ns-a))


# all_locations have the same lambda --------------------------------------
eta = matrix(0,nrow = ns, ncol = nt)
eta[idx1, ] = matrix(beta_1 %*% t(tPHI))
eta[idx2, ] = matrix(beta_2 %*% t(tPHI))
matplot(t(eta), xlab = 't', ylab = 'eta', main = 'Simulated linear predictor')

mu = matrix(0,nrow = ns, ncol = nt)
for(i in 1:ns){mu[i,] = exp(eta)[i,]*boundary$population[i]}

Y = rpois(length(mu),lambda = mu)

Y_obs = matrix(Y, nrow = ns, ncol = nt)
ratio = matrix(0,nrow = ns, ncol = nt)
for(i in 1:ns){ratio[i,] = as.matrix(Y_obs[i,]/boundary$population[i])}
raw_data = cbind(coords, cluster_true,boundary$population,Y_obs)
# factor = c(1:27)[idx1]
 matplot(t(Y_obs/boundary$population), type = "b",pch=1,xlab = 't', ylab = 'Y/pop', main = 'Observed linear predictor')


save(raw_data,file = 'data/state_data_sim2_eqlambda_eqvar.Rdata')
save(boundary,file = '~/Documents/project3/data/boundary.Rdata')

# ggplot(data = boundary) +
#   geom_sf()


# Functions ---------------------------------------------------------------

plotTrueMean<- function(beta_true, tPHI){
  eta_true = beta_true %*% t(tPHI)
  matplot(exp(t(eta_true)), type = "b",pch=1,col = 1:27)
}
