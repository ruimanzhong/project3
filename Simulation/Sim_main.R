#' @  k_true number of cluster 
#' @ para ns number of locations 
#' @ para nt number of time points
library(RhpcBLASctl)

blas_get_num_procs()
blas_set_num_threads(6)
rm(list = ls())
packages = c("wavethresh", "gridExtra","tidyverse",
             "RColorBrewer",'INLA', "deldir", "fields", "igraph",
             "ggraph", "RANN" , "grid", "parallel", "class", "sf", 
             "ggplot2",'eSDM', 'fda')
## Now load or install&load all packages
package.check <- lapply( packages,
                         FUN = function(x) {
                           if (!require(x, character.only = TRUE)) { install.packages(x, dependencies = TRUE) 
                             library(x, character.only = TRUE)
                           }else
                             library(x, character.only = TRUE)
                         } )
library(Matrix.utils)
 # source('~/Documents/BSTDM/project3/Simulation/Sim_functions.R')
source('~/Documents/BSTDM/project3/Simulation/Sim2_function.R')
# Sim data loadin--------------------------------------------------------------------
# Sim 1
load("~/Documents/project3/data/state_data_novar_lambda1.Rdata")
# Sim 2
load("~/Documents/BSTDM/project3/Simulation/data/state_data_sim2_eqlambda_eqvar.Rdata")
# Built initial tree ------------------------------------------------------
coords = cbind(raw_data[,1], raw_data[,2])
colnames(coords) = c("lon", "lat")
cluster_true= raw_data[,3]
population = raw_data[,4]
Y = raw_data[,5:132]
Y_obs = Y
ns = dim(Y)[1]
nt = dim(Y)[2] 
ratio = matrix(0,nrow = ns, ncol = nt)
for(i in 1:ns){ratio[i,] = as.matrix(Y_obs[i,]/population[i])}
t = seq(0,1, length.out = nt)
J = 3
##initial graph
library(rgeoboundaries)
boundary <- geoboundaries("Brazil","adm1")
ggplot(data = boundary) +
  geom_sf()

adjacency_matrix <- st_touches(boundary) %>% as("matrix")

# Convert the matrix to igraph graph object
graph <- graph.adjacency(adjacency_matrix, mode = "undirected", diag = FALSE)
graph0 = graph

###Inverse discrete wavelet transform
# splinebasis = create.bspline.basis(c(0,1), J)
# tPHI= eval.basis(t, splinebasis) 

# polynomial basis --------------------------------------------------------
J = 3
tPHI = matrix(c(rep(1,nt), t, t^2, t^3), nrow = nt)
dim(tPHI)

 # iteration procedure -----------------------------------------------------
##Initial beta estimates
beta.ini= ginv(t(tPHI) %*% tPHI) %*% t(tPHI) %*% t(log(ratio)) # p * s
beta.ini = t(beta.ini) 

inc.mat=get.edgelist(graph0, names = F)
weights=sqrt(apply((beta.ini[inc.mat[,1],]-beta.ini[inc.mat[,2],])^2,1,sum))
E(graph0)$weights=weights
mstgraph.ini=mst(graph0)


V(mstgraph.ini)$vid=1:ns

##Create 7 initial clusters
##id of vertices for removal
 # rmid=order(E(mstgraph.ini)$weights, decreasing = TRUE)[1:2]
graph_comp=components(delete.edges(mstgraph.ini, c(1,4)))
membership.ini=graph_comp$membership
membership.ini
membership=membership.ini ##cluster membership
mstgraph=mstgraph.ini     ## initial MST


# Initial parameters

 Y=y;  ##Specify the response matrix
 # X = ginv(tPHI) # 
# J=log2(nt) ## The finest resolution of the wavelet representation 

#####Initial values of beta and sigma2
clust_uniq=sort(unique(membership.ini)); k=length(clust_uniq)
beta_uniq.ini=matrix(0, k, ncol(tPHI))

sigmasq.ini=rep(0.01, k)

# for(i in 1:k){
#   ns1=length(which(membership.ini==i))
#   
#   log_y_clust=matrix(log(y)[which(membership.ini==i),], ns1, nt)
#   if(ns1>1){
#   beta_uniq.ini[i,]=t(as.matrix(colSums(1/ns1*beta.ini[which(membership.ini==i),])))
#   } else {
#   beta_uniq.ini[i,]=t(as.matrix(beta.ini[which(membership.ini==i),]))
#   }
#   #Heterogeneous variances
#   if(ns1>2){
#     sigmasq.ini[i]=sum((log_y_clust-kronecker(matrix(1, ns1, 1),
#                                           beta_uniq.ini[i,]%*%t(tPHI)))^2)/(ns1*nt-ncol(tPHI))
#   }else{
#     ##If a cluster only has one observation, give an initial value
#     sigmasq.ini[i]=0.01
#   }
#   
#   
# }
###Set the initial values of lambda, pi, and indicator variables gamma.
lambda.ini=matrix(1, k, J+1)
# para_pi.ini=matrix(0.5, k, J+1)
# para_pi.ini[, 1]=1
# ind_gam.ini=matrix(1, k, nt)
X = matrix(c(t, t^2, t^3), nrow = nt)
tPHI = X
###List of initial values
init_val<-list("beta"=beta_uniq.ini, "cluster"=membership.ini, 
               "sigmasq_y"=sigmasq.ini, 
               "lambda"=lambda.ini, 
               # "para_pi"=para_pi.ini, "ind_gam"=ind_gam.ini, 
               "trees"=mstgraph.ini)

hyperpar<-list("a0"= 0.01, "b0"= 0.01, "c0"=NULL, "d0"= NULL, "c"= 0.9)

##THIN: thinning by recording every (THIN)-th iteration
MCMC=2000; BURNIN=50; THIN=5
n.cores=6

###path and file name for saving the results
path_save="~/Documents/BSTDM/project3/Simulation/results/Sim2_MCMC_demo.RData"

start_time <- Sys.time()
Fun_Wave_Clust(Y, tPHI, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, path_save, seed = 1234, n.cores = 6)
end_time <- Sys.time()

clust_idx = 100
clust_max=cluster_out[clust_idx,]
mstgraph=MST_out[[clust_idx]]


title1="true clusters"
title2="Fclust-RST"
colnames(coords)=c("lon", "lat")
p1=plotGraphData(coords, mstgraph, cluster_true, title1)
p2=plotGraphData(coords, mstgraph, clust_max, title2)
grid.arrange(p1, p2, nrow = 1)

Y_obs[which(clust_max == 3),]/population[which(clust_max == 3)]
ratio = matrix(0,nrow = ns, ncol = nt)
for(i in 1:ns){ratio[i,] = as.matrix(Y_obs[i,]/population[i])}
plot(t,Y_obs[which(clust_max == 3),]/population[which(clust_max == 3)])
lines(t,ratio[15,], lch = 3)

matplot(t(ratio* 100000))


matplot(t(ratio[which(cluster_true == 1),]), type = 'l')
ratio
