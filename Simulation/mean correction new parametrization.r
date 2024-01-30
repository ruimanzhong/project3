library(parallel)
library(Matrix)
library(statmod)
library(runjags)
library(rjags)
library(INLA)
library(sn)
library(mvtnorm)

set.seed(12345)
n <- 50
beta <- -1.75
u <-  rep(beta, n) + rep(-0.386832, n)   #rnorm(n, mean = 0, sd = 1) 
mu <- 1
y <- rpois(n, mu)
plot(1:n, y, ylim = c(min(min(y), min(mu)), max(max(y), max(mu))), main = "Poisson Observation")
lines(1:n, mu, col = "red")
#lines(1:n, x.star[1:n], col = "blue")


#INLA Result
idx = 1:n
prior.ui <- list(prec = list(prior = "gaussian", param = c(0, 1)))
formula =y ~ 1+f(idx,model="iid",hyper = prior.ui)

res = inla(formula = formula,
           family = "poisson",
           data = list(y=y, idx = idx),
           control.compute = list(config=T, return.marginals.predictor=T),
           control.predictor = list(compute = T),
           control.fixed=list(prec.intercept = 1),
           control.inla = list(control.vb=list(enable = TRUE)),
           inla.mode = "experimental")
#res = inla.hyperpar(res)

# log of the distribution of the log precision
# lprior.theta <- function(theta){
#   return(theta - exp(theta))
# }

# # log of the distribution of the log precision
lprior.theta <- function(theta){
  return(dnorm(theta, mean = 0, sd = 1, log = T))
}

# A
A <-  cbind(diag(1, n), rep(1, n))
  
  
# GMRF
Q <- function(theta){
  #Q <- matrix(0, nrow = 2*n, ncol = 2*n)
  Q <- diag(c(rep(exp(theta), n), 1))
  Q <- as(Q, "sparseMatrix")
  return(Q)
}

# Prior of X given theta-

def.xprior <- function(theta){
  return(list(Q = Q(theta), b = rep(0, n + 1)))
}

# log of Normal density
log.normal <- function(x, def){
  def <- def.xprior(1.23)
  mu <- solve(def$Q, def$b)
  x <- x - mu
  q <- -0.5*crossprod(x ,def$Q)%*%x
  return(drop(q) + sum(log(diag(chol(def$Q)))))
}

# distribution of x given y and theta
def.xcond <- function(theta){
  initial <- c(log(1+y), 0)
  for (i in 1:10) {
    initial.1 <- A%*%initial
    initial.1 <- as.numeric(initial.1)
    D <- diag(exp(initial.1))
    H <- def.xprior(theta)$Q + crossprod(A, D)%*%A
    b <- y - exp(initial.1) + initial.1*exp(initial.1)
    b <- crossprod(b, A)
    updated <- solve(H, as.numeric(b))
    
    
    #if (T) print(mean(abs(initial - updated)))
    if(mean(abs(initial - updated)) < 1e-8)
      break
    initial <- updated
  }
   x <- as.numeric(updated) 
   x.cond <- list(Q = H, b = t(b), x = x)
  # return(list(H = H, b = b, x = x))
  #  stop("")
  return(x.cond)
}
# Log Posterior of theta
log.post <- function(theta){
  initial <- c(log(1+y), 0)
  prior <- def.xprior(theta)
  for (i in 1:10) {
    initial.1 <- A%*%initial
    initial.1 <- as.numeric(initial.1)
    D <- diag(exp(initial.1))
    H <- prior$Q + crossprod(A, crossprod(D, A))
    b <- y - exp(initial.1) + initial.1*exp(initial.1)
    b <- crossprod(b, A)
    updated <- solve(H, as.numeric(b))
    
    
    #if (T) print(mean(abs(initial - updated)))
    if(mean(abs(initial - updated)) < 1e-8)
      break
    initial <- updated
  }
  x <- as.numeric(updated) 
  x.cond <- list(Q = H, b = b)
  # return(list(H = H, b = b, x = x))
  #  stop("")
  
  D <- as.numeric(A%*%x)
  ldens <- lprior.theta(theta) + log.normal(x, def.xprior(theta)) + 
    sum(dpois(y, lambda = exp(D), log = TRUE)) - sum(log(diag(chol(x.cond$Q))))
  return(drop(ldens))   
}

grid.interval <- 0.1
theta.s <- seq(-3, 3, by = grid.interval)
#theta.s <- unlist(lapply(1:nconfigs,FUN = function(i){res$misc$configs$config[[i]]$theta}))
ldens <- unlist(mclapply(theta.s, log.post, mc.cores = detectCores()))
ldens <- exp(ldens - max(ldens))

ldens <- ldens/(sum(ldens)*grid.interval)
plot(theta.s, ldens, type = "l")

# npts <- res$misc$configs$nconfig
# theta.s <- as.numeric(unlist(mclapply(1:npts, function(i) return(res$misc$configs$config[[i]]$theta), mc.cores = detectCores())))
# #theta.s <- unlist(lapply(1:nconfigs,FUN = function(i){res$misc$configs$config[[i]]$theta}))
# ldens <- unlist(mclapply(theta.s, log.post, mc.cores = detectCores()))
# ldens.store <- ldens
# ldens <- exp(ldens - max(ldens))
# 
# ldens <- ldens/(sum(ldens[1:(npts-1)]*diff(theta.s)))
# plot(theta.s, ldens, type = "l")

#JAGS
model = "model {
  for(i in 1:n){
    y[i] ~ dpois(lambda[i])
    log(lambda[i]) = beta + u[i]
    u[i] ~ dnorm(0, tau)
  }

  beta ~ dnorm(0, 1)
  tau ~ dnorm(0, 1)
  theta <- log(tau)
}"

res.1 = run.jags(model = model,
                 monitor = c("tau", "theta"),
                 data = list(y = y, n = n),
                 n.chains = 1,
                 inits = list(tau = 1),
                 burnin = 1000,
                 sample = 25000)

trace = combine.mcmc(res.1)
hist(trace[,2], n = 200, probability = T, col = "gray")
lines(theta.s, ldens, type = "l", col = "red", lwd = 1)
post.prec <- res$marginals.hyperpar$`Precision for idx`[,1]
post.prec.dens <-res$marginals.hyperpar$`Precision for idx`[,2]
post.log.prec <-log(post.prec)
post.log.prec.dens  = post.prec.dens*post.prec
lines(post.log.prec, post.log.prec.dens, lwd = 1, col = "blue")

#Smoothing the curve of posterior distribution
theta.dist <- splinefun(theta.s, ldens)
hist(trace[,2], n = 200, probability = T, col = "gray")
lines(seq(-2, 2.5, by = 0.001), theta.dist(seq(-2, 2.5, by = 0.001)), type = "l", col = "red", lwd = 1)
lines(post.log.prec, post.log.prec.dens, lwd = 1, col = "blue")

# Normal Approximation to the posterior
# r.post.theta <- optim(0, theta.dist, NULL, control = list(fnscale = -1), method = "BFGS")
# theta.star <- r.post.theta$par
# 
# theta.dist.1 <- function(theta) {
#   if (theta.dist(theta) < 0) {
#     return(1e-20)
#   }else return(theta.dist(theta))
# }
# 
# KLD <- function(delta){
#   out <- gauss.quad.prob(21, dist = "normal", mu = delta[1], sigma = exp(delta[2]))
#   store <- unlist(mclapply(out$nodes, theta.dist.1, mc.cores = detectCores()))
#   store <- log(store)
#   return(sum(out$weights*(dnorm(out$nodes, mean = delta[1], sd = exp(delta[2]), log = T) -
#                             store)))
# }
# 
# r.KLD <- optim(c(0,0), KLD, NULL, method = "BFGS")
# 
# hyp.mean <- r.KLD$par[1]
# hyp.sd <- exp(r.KLD$par[2])
# 
# plot(seq(-2, 2.5, 0.01), theta.dist(seq(-2, 2.5, 0.01)), type = "l", lwd = 2, col = "red", main = "Distribution of theta",  xlab = expression(theta), ylab = expression(f(theta)))
# lines(seq(-2, 2.5, 0.01), dnorm(seq(-2, 2.5, 0.01),mean = hyp.mean, sd = hyp.sd), type = "l", lwd = 2, col = "darkgreen", lty = 2)
# lines(seq(-2, 2.5, 0.01), exp(lprior.theta(seq(-2, 2.5, 0.01))), col = "blue", lwd = 2)
# 
# legend(x = "topright",
#        legend = c("prior", "INLA", "N.Approx Post"),
#        lty = c(1, 1, 2),
#        col = c("blue", "red", "darkgreen"),
#        lwd = 2)


# x.star <- solve(def.xcond(theta.star)$Q, t(def.xcond(theta.star)$b))
# plot(1:n, y, ylim = c(min(min(y), min(mu)), max(max(y), max(mu))))
# lines(1:n, mu, col = "red")
# plot(1:n, x.star[1:n], col = "blue", type = "l")
# plot(1:n, res$summary.linear.predictor$mean, col = "green")

# poillike <- function(y, mu, sigsq){
#   like <- exp(mu+(0.5*sigsq)) - y*mu
#   return(like)
# }
# 
# 
# xcond.c <- function(theta){
#   store1 <- def.xcond(theta)
#   store0 <- def.xprior(theta)
# 
#   VB <- function(delta){
#     b.c <- store1$b
#     b.c[n+1] <- b.c[n+1] + delta
#     mu1 <- solve(store1$Q, b.c)
#     SIGMA1 <- solve(store1$Q)
#     KLD <- 0.5*(sum(diag(crossprod(store0$Q, SIGMA1)))
#                 - (n+1)
#                 + crossprod(mu1, crossprod(store0$Q, mu1))
#                 + sum(log(diag(chol(store1$Q)))) - sum(log(diag(chol(store0$Q)))))
#     # poillike <- function(y, mu, sigsq){
#     #   like <- exp(mu+(0.5*sigsq)) - y*mu
#     #   return(like)
#     # }
#     llike <- unlist(mcmapply(function(a1, a2, a3) return(poillike(a1, a2, a3)), a1 = y, a2 = A%*%mu1, a3 = diag(A%*% tcrossprod(SIGMA1, A)),
#                              mc.cores = detectCores()))
#     llike <- sum(llike)
#     return(llike + drop(KLD))
#   }
# 
#   r <- optim(0, VB, NULL, method = "BFGS")
#   bcorrc <- store1$b
#   bcorrc[n+1] <- bcorrc[n+1] + r$par
#   x <- solve(store1$Q, bcorrc)
#   return(list(Q = store1$Q, b = bcorrc, x = x))
# }
# 
# 
# THETA <- function(theta){
#   store <- def.xcond(theta)
#   x.theta <- A%*%as.numeric(store$x)
#   store1 <- (sum(log(diag(chol(Q(theta)))))
#              + sum((-rep(1, n) + x.theta - 0.5*x.theta^2)*exp(x.theta))
#              + 0.5*crossprod(store$b, solve(store$Q,store$b))
#              - sum(log(diag(chol(store$Q)))))
#   return(drop(-store1))
# }
# 
# 
# eval.point <- 15
# VB.1 <- function(delta){
#   #delta <- c(0,0)
#   out <- gauss.quad.prob(eval.point, dist = "normal", mu = hyp.mean + delta[1],
#                          sigma = hyp.sd*exp(delta[2]))
#   KLD <- sum(out$weights*(dnorm(out$nodes, mean = hyp.mean + delta[1],
#                                 sd = hyp.sd*exp(delta[2]), log = T) -
#                             lprior.theta(out$nodes)))
#   node <- unlist(mclapply(out$nodes, THETA, mc.cores = detectCores()))
#   node <- as.numeric(node)
#   llike <- sum(out$weights*node)
#   #llike <- 0
#   return(llike + KLD)
# }
# 
# r.1 <- optim(c(0, 0), VB.1, NULL, method = "BFGS")
# 
# hist(trace[,2], n = 200, probability = T, col = "gray", main = "Distribution of theta",
#      xlab = expression(theta), ylab = expression(f(theta)))
# lines(seq(-2, 2.5, 0.01), theta.dist(seq(-2, 2.5, 0.01)), type = "l", lwd = 2, col = "red")
# #lines(seq(-2, 2.5, 0.01), dnorm(seq(-2, 2.5, 0.01), mean = hyp.mean, sd = hyp.sd), type = "l", lwd = 2, col = "red", lty = 1)
# lines(seq(-2, 2.5, 0.01), dnorm(seq(-2, 2.5, 0.01), mean = hyp.mean + r.1$par[1], sd = hyp.sd*exp(r.1$par[2])), type = "l", lwd = 2, col = "darkgreen", lty = 2)
# lines(seq(-2, 2.5, 0.01), exp(lprior.theta(seq(-2, 2.5, 0.01))), col = "deeppink", lwd = 2, lty = 3)
# legend(x = "topright",
#        legend = c("prior", "INLA", "VB"),
#        lty = c(3, 1, 1),
#        col = c("deeppink", "red", "darkgreen"),
#        lwd = 2)
# 

# log.post.c <- function(theta){
#   store1 <- def.xcond(theta)
#   store0 <- def.xprior(theta)
# 
#   VB <- function(delta){
#     b.c <- store1$b
#     b.c[(n+1)] <- b.c[(n+1)] + delta
#     mu1 <- solve(store1$Q, b.c)
#     SIGMA1 <- solve(store1$Q)
#     KLD <- 0.5*(sum(diag(store0$Q%*%SIGMA1))
#                 - (n+1)
#                 + crossprod(mu1, store0$Q)%*%mu1
#                 + sum(log(diag(chol(store1$Q)))) - sum(log(diag(chol(store0$Q)))))
#     # poillike <- function(y, mu, sigsq){
#     #   like <- exp(mu+(0.5*sigsq)) - y*mu
#     #   return(like)
#     # }
#     llike <- unlist(mcmapply(function(a1, a2, a3) return(poillike(a1, a2, a3)), a1 = y, a2 = A%*%mu1, a3 = diag(A%*% tcrossprod(SIGMA1, A)),
#                              mc.cores = detectCores()))
#     llike <- sum(llike)
#     return(llike + drop(KLD))
#   }
# 
#   r <- optim(0, VB, NULL, method = "BFGS")
#   bcorrc <- store1$b
#   bcorrc[(n+1)] <- bcorrc[(n+1)] + r$par
#   x <- solve(store1$Q, bcorrc)
#   ldens <- lprior.theta(theta) + log.normal(x, def.xprior(theta)) +
#     sum(dpois(y, lambda = exp(x[1:n]), log = TRUE)) - sum(log(diag(chol(store1$Q))))
#   return(drop(ldens))
# }




# Log Posterior of theta egil
ldens.c <- numeric(npts)
for (i in 1:npts) {
  impv <- as.matrix(res$misc$configs$config[[i]]$mean - res$misc$configs$config[[i]]$improved.mean)
  prec.mat <- res$misc$configs$config[[i]]$Q
  prec.mat <- prec.mat + t(prec.mat)
  diag(prec.mat) <- diag(prec.mat)/2
  ldens.c[i] <- ldens.store[i] + drop(0.5*crossprod(impv, crossprod(prec.mat, impv)))
}

ldens.c[1] <- -42
ldens.c <- exp(ldens.c - max(ldens.c))

ldens.c <- ldens.c/(sum(ldens.c[1:(npts-1)]*diff(theta.s)))
plot(theta.s, ldens.c, type = "l")


# grid.interval <- 0.5
# theta.s <- seq(-3, 3, by = grid.interval)
# #theta.s <- unlist(lapply(1:nconfigs,FUN = function(i){res$misc$configs$config[[i]]$theta}))
# ldens <- unlist(mclapply(theta.s, log.post, mc.cores = detectCores()))
# ldens <- exp(ldens - max(ldens))
# 
# ldens <- ldens/(sum(ldens)*grid.interval)
# plot(theta.s, ldens, type = "l")



hist(trace[,2], n = 200, probability = T, col = "gray")
lines(theta.s, ldens, type = "l", col = "red", lwd = 2, lty = 1)
lines(theta.s, ldens.c, type = "l", col = "darkgreen", lwd = 2, lty = 2)
legend("topright",
       legend = c("INLA", "VB"),
       lty = c(1, 2),
       lwd = c(2,2),
       col = c("red", "darkgreen"))


# ldens.c <- numeric(npts)
# for (i in 1:npts) {
#   theta <- as.numeric(res$misc$configs$config[[i]]$theta)
#   #x <- res$misc$configs$config[[i]]$mean
#   x <- (res$misc$configs$config[[i]]$improved.mean)
#   precison.mat <- res$misc$configs$config[[i]]$Q
#   precison.mat <- precison.mat + t(precison.mat)
#   diag(precison.mat) <- diag(precison.mat)/2
#   ldens.c[i] <- (lprior.theta(theta) + log.normal(x, def.xprior(theta)) +
#     sum(dpois(y, lambda = exp(x[1:n]), log = TRUE)) - sum(log(diag(chol(precison.mat)))))
# }
# 
# ldens.c <- exp(ldens.c - max(ldens.c))
# ldens.c <- ldens.c/(sum(ldens.c[1:(npts-1)]*diff(theta.s)))
# plot(theta.s, ldens.c, type = "l")
# 
# hist(trace[,2], n = 200, probability = T, col = "gray")
# lines(seq(-2, 2.5, by = 0.001), theta.dist(seq(-2, 2.5, by = 0.001)), type = "l", col = "red", lwd = 1)
# lines(theta.s, ldens.c, type = "l")

# tik <- Sys.time()
# ldens.c <- unlist(mclapply(theta.s, log.post.c, mc.cores = detectCores()))
# tok <- Sys.time()
# print(tok - tik)
# #DELTA <- ldens.c[2*(1:length(theta.s))]
# #ldens.c <- ldens.c[2*(1:length(theta.s))-1]
# ldens.c <- exp(ldens.c - max(ldens.c))
# ldens.c <- ldens.c/(sum(ldens.c)*grid.interval)
# plot(theta.s, ldens.c, type = "l")
# #plot(theta.s, DELTA, type = "l")
# 
# hist(trace[,2], n = 200, probability = T, col = "gray", ylim = c(0,2))
# lines(theta.s, ldens, type = "l", col = "red", lwd = 2, lty = 1)
# lines(theta.s, ldens.c, type = "l", col = "darkgreen", lwd = 2, lty = 2)
# legend("topright",
#        legend = c("INLA", "VB"),
#        lty = c(1, 2),
#        lwd = c(2,2),
#        col = c("red", "darkgreen"))
# 
# #save.image("result.RData")
