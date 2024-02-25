library("splines")
library(INLA)
knots = seq(0,1, length.out  = J)
n <- 128
beta_0 <- -5
beta_1 <- 3
t <- seq(0,1, length.out = n)
population = 10000

u1 <-  exp(rnorm(n, mean = beta_0, sd = 0.1) + rnorm(n, mean = beta_1, sd = 1)*t ) *population
set.seed(12345)
# u2 = exp(rnorm(n, mean = beta_0, sd = 0.1) + rnorm(n, mean = beta_1, sd = 1)*t) *population
y1 <- rpois(n, u1)
set.seed(123)
y2 <- rpois(n, u1)
y <- c(y1, y2)
plot(1:n, y1/population, ylim = c(min(min(u1/population), min(mu)), max(max(u1/population), max(mu))), type = 'l', main = "Poisson Observation")
#lines(1:n, mu, col = "red")
lines(1:n, y2/population, col = "blue")

idx = rep(1:n,2)
p = 2
prior.ui <- list(prec = list(prior = "loggamma", param = c(0.01, 0.01)))
formula =y ~ 1+t+f(idx,model="iid",hyper = prior.ui)

res = inla(formula = formula,
           family = "poisson",E = population,
           data = list(y=y, t = rep(t,2), idx = idx),
           control.compute = list(config=T, return.marginals.predictor=T,mlik = T),
           control.predictor = list(compute = T),
           control.fixed=list(prec.intercept = 1,prec = 1),
           control.inla = list(control.vb=list(enable = TRUE)),
           inla.mode = "experimental")
summary(res)
res1 = inla(formula = formula,
           family = "poisson",E = population,
           data = list(y=y, t =rep(t,2), idx = 1:256),
           control.compute = list(config=T, return.marginals.predictor=T,mlik = T),
           control.predictor = list(compute = T),
           control.fixed=list(prec.intercept = 1, prec = 1),
           control.inla = list(control.vb=list(enable = TRUE)),
           inla.mode = "experimental")

summary(res1)

res2 = inla(formula = formula,
            family = "poisson",E = population,
            data = list(y=y2, t =t, idx = idx[1:128]),
            control.compute = list(config=T, return.marginals.predictor=T,mlik = T),
            control.predictor = list(compute = T),
            control.fixed=list(prec.intercept = 1, prec = 1),
            control.inla = list(control.vb=list(enable = TRUE)))

summary(res2)
 res$mlik[[1]] -(res1$mlik[[1]] + res1$mlik[[2]])
# # log of the distribution of the log precision

lprior.theta <- function(theta){
  return(dgamma(theta, shape = 0.01, scale = 0.01, log = T))
}


# GMRF
Q <- function(theta){
  #Q <- matrix(0, nrow = 2*n, ncol = 2*n)
  Q <- diag(c(exp(theta),rep(1,p)))
  Q <- as(Q, "sparseMatrix")
  return(Q)
}

# Prior of X given theta-

def.xprior <- function(theta){
  return(list(Q = Q(theta), b = rep(0, n + p)))
}

# distribution of x given y and theta
def.xcond <- function(A, y, theta, population, n){
  initial <- c(log(1+y[1:n])/population,0, 0)
  for (i in 1:100) {
    initial.eta <- A%*%initial
    initial.eta <- as.numeric(initial.eta)
    D <- diag(exp(initial.eta)*population)
    H <- def.xprior(theta)$Q + crossprod(A, D)%*%A
    b <- y - population*exp(initial.eta) +population*initial.eta*exp(initial.eta) 
    b <- crossprod(b,A)
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
# A
ns = 2
A.single <-  cbind(diag(1, n), t, rep(1, n))
A = do.call(rbind, replicate(ns, A.single, simplify = FALSE))
theta <- rgamma(128,shape = 0.01, scale = 0.01)
x.prior = def.xprior(theta)
x.cond = def.xcond(A, y, theta, population, n)
# conditional marginal likelihood -----------------------------------------f
mlik.cond <- function(A, y, x.cond, x.prior){
 x.mode <- x.cond$x
 eta.mode <- A %*% x.mode
 log.x.prior <-   - 0.5*crossprod(x.mode, x.prior$Q)%*%x.mode 
 log_likelihood = crossprod(y, eta.mode) + sum(exp(eta.mode))
 log_x.cond = - 0.5*crossprod(x.mode, x.cond$Q)%*%x.mode + crossprod(x.cond$b,x.mode) 
 mlik <- log.x.prior + log_likelihood - log_x.cond
 return(mlik)
}                                            
mlik <- mlik.cond(A, y, x.cond, x.prior)
