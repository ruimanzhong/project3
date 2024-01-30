#' @  k_true number of cluster 
#' @ para ns number of locations 
#' @ para nt number of time points


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
source('~/Documents/project3/Code_Fclust-RST/SCF_wave_functions.R')
load("~/Documents/project3/data/raw_data.Rdata")
dim(raw_data)
# Built initial tree ------------------------------------------------------
coords = cbind(raw_data[,1], raw_data[,2])
colnames(coords) = c("lon", "lat")
cluster_true = raw_data[,3]
population = raw_data[,4]
y = raw_data[,5:260]/population
ns = dim(y)[1]
nt = dim(y)[2] 
t = seq(0,1, length.out = nt)
J = 8
##initial graph
library(rgeoboundaries)
boundary <- geoboundaries("Brazil","adm2")
ggplot(data = boundary) +
  geom_sf()

adjacency_matrix <- st_touches(boundary) %>% as("matrix")

# Convert the matrix to igraph graph object
graph <- graph.adjacency(adjacency_matrix, mode = "undirected", diag = FALSE)
graph0 = graph

###Inverse discrete wavelet transform
splinebasis = create.bspline.basis(c(0,1), J)
tPHI= eval.basis(t, splinebasis) 



# iteration procedure -----------------------------------------------------
##Initial beta estimates
beta.ini= ginv(t(tPHI) %*% tPHI) %*% t(tPHI) %*% t(log(y)) # p * s
beta.ini = t(beta.ini) 

inc.mat=get.edgelist(graph0, names = F)
weights=sqrt(apply((beta.ini[inc.mat[,1],]-beta.ini[inc.mat[,2],])^2,1,sum))
E(graph0)$weights=weights
mstgraph.ini=mst(graph0)

V(mstgraph.ini)$vid=1:ns

##Create 7 initial clusters
##id of vertices for removal
rmid=order(E(mstgraph.ini)$weights, decreasing = TRUE)[1:8]

graph_comp=components(delete.edges(mstgraph.ini, rmid))
membership.ini=graph_comp$membership

membership=membership.ini ##cluster membership
mstgraph=mstgraph.ini     ## initial MST

plotGraph(coords, mstgraph.ini, 'Initial MST') + geom_point(aes(x = lon, y = lat, color = factor(cluster_true)),
                                                        size = 1, data = as.data.frame(coords))
# Initial parameters

 Y=y;  ##Specify the response matrix
 # X = ginv(tPHI) # 
# J=log2(nt) ## The finest resolution of the wavelet representation 

#####Initial values of beta and sigma2
clust_uniq=sort(unique(membership.ini)); k=length(clust_uniq)
beta_uniq.ini=matrix(0, k, ncol(tPHI))

sigmasq.ini=rep(0.01, k)

for(i in 1:k){
  ns1=length(which(membership.ini==i))
  
  log_y_clust=matrix(log(y)[which(membership.ini==i),], ns1, nt)
  if(ns1>1){
  beta_uniq.ini[i,]=t(as.matrix(colSums(1/ns1*beta.ini[which(membership.ini==i),])))
  } else {
  beta_uniq.ini[i,]=t(as.matrix(beta.ini[which(membership.ini==i),]))
  }
  #Heterogeneous variances
  if(ns1>2){
    sigmasq.ini[i]=sum((log_y_clust-kronecker(matrix(1, ns1, 1),
                                          beta_uniq.ini[i,]%*%t(tPHI)))^2)/(ns1*nt-ncol(tPHI))
  }else{
    ##If a cluster only has one observation, give an initial value
    sigmasq.ini[i]=0.1
  }
  
  
}
###Set the initial values of lambda, pi, and indicator variables gamma.
lambda.ini=matrix(1, k, J+1)
# para_pi.ini=matrix(0.5, k, J+1)
# para_pi.ini[, 1]=1
# ind_gam.ini=matrix(1, k, nt)

###List of initial values
init_val<-list("beta"=beta_uniq.ini, "cluster"=membership.ini, 
               "sigmasq_y"=sigmasq.ini, 
               "lambda"=lambda.ini, 
               # "para_pi"=para_pi.ini, "ind_gam"=ind_gam.ini, 
               "trees"=mstgraph.ini)

hyperpar<-list("a0"=4, "b0"=2*mean(sigmasq.ini), "c0"=rep(4, J+1), "d0"=rep(2, J+1), "c"= 0.5)

##THIN: thinning by recording every (THIN)-th iteration
MCMC=20000; BURNIN=5000; THIN=5
n.cores=2

###path and file name for saving the results
path_save="simulation/Sim_MCMC_demo.RData"

Fun_Wave_Clust(log(Y), X, graph0, init_val, hyperpar, MCMC, BURNIN, THIN, path_save, seed = 1234, n.cores)
