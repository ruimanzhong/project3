packages = c( "gridExtra","tidyverse",
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

# J is the number of basis function, J+1 includes intercept
# k is the number of clusters
# X is the design matrix without intercept

k = 3
J = 3
X = matrix(c(t, t^2, t^3), nrow = nt)

sigmasq_y = rep(0.01, k)
lambda=matrix(1, k, J+1) 

# Y is a ns times nt matrix
# ind control the memberships to fit the model
# population is the population data

mlik <- function(ind){
  ns=length(ind)
  Yk= matrix(Y[ind,],ncol = nt)
  vec_Yk = as.integer(as.vector(t(Yk))) 
  
  sigma2_yk=sigmasq_y[k]
  
  if(p>1){
    lambda_k=lambda[k,]
  }else{
    lambda_k=lambda
  }
  COV = cbind(1, do.call(rbind, replicate(nrow(Yk), X, simplify = FALSE)))
  beta.prec = diag(1)*(1/sigma2_yk)
  data = as.data.frame(cbind(vec_Yk = vec_Yk,COV = COV))
  
  idx = rep(1:nt,ns)
  formula= vec_Yk~  0 + COV + f(idx, model = 'iid',
                            hyper = list(prec = list(prior = "loggamma", param = c(a0,b0))))
  m.bs3 <- INLA::inla(formula, family = "poisson", E = rep(population[ind],each = nt),
                      data = data, 
                      control.fixed = list(prec = beta.prec, prec.intercept = 0.38)
  )
  print(summary(m.bs3))
  return(m.bs3$mlik[[1]] )
}
ind_all = which(cluster_true == 1)
# merge step

 ind_2 = ind_all[c(1:7,9:11)]
ind_1 = ind_all[8]
M = mlik(ind_all)
M1 = mlik(ind_1)
M2 = mlik(ind_2)
 M-(M1+M2)
 
