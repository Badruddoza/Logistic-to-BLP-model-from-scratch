rm(list=ls())
#BLP 1995 replication
#Following Paul Schrimpf 
#uses mean price not log price
if (!require(hdm)) install.packages("hdm")
library(hdm)
data(BLP)
tab2 <- c(3.393, 6.711, 8.728, 13.074, 68.597)
stopifnot(all(abs(quantile(BLP$price, c(0,0.25,0.5,0.75,1)) - tab2 +
                    11.761)<0.005))
BLP$price <- BLP$price + 11.761

share.fn <- function(delta,  ## J vector
                     x,      ## J by K
                     log.y,  ## S vector of log(y_i)
                     log.yp, ## S by J of log(y_i - p_j) 
                     v,      ## S by K
                     alpha,  ## scalar
                     sigma)  ## K vector
{
  stop("## TODO: compute J vector of shares in this market")
  return(share)
}

delta.fn <- function(s,x,log.y,log.yp,v,alpha,sigma,
                     tol=1e-6, tol.s=1e-8,
                     max.iter=100)
{
  delta.new <- log(s) - log(1-sum(s)) ## initial guess
  dd <- 1 ## ||change in delta||
  ds <- 1 ## ||observed shares - share.fn(delta)||
  iter <- 0
  while ((dd>tol) | (ds>tol.s)) {
    delta.old <- delta.new
    delta.new <- stop("## TODO: update delta using contraction mapping")
    dd <- max(abs(delta.old - delta.new))
    ds <- max(abs(s - sm))
    iter = iter+1
    if (iter>max.iter) {
      warning(sprintf("Maximum iterations (%d) reached, returning with norm(delta.new - delta.old) = %.2g", 
                      max.iter, dd))
      break;
    }
  }
  return(delta.new)
}

## Create instruments
## demand instruments
Z <- hdm:::constructIV(BLP$firm.id, BLP$cdid, BLP$id,
                       cbind(1,BLP[,c("hpwt","air","mpd","mpg","space","price")]))
## supply instruments
W <- log(BLP[,c("hpwt","mpg","space","mpd")])
colnames(W) <- paste("log",colnames(W),sep=".")
Wiv <- hdm:::constructIV(BLP$firm.id, BLP$cdid, BLP$id, W)

## Draws for simluting integral
S <- 100 
T <- length(unique(BLP$cdid))
K <- 5
set.seed(41658)
y.s <- matrix(exp(rnorm(S*T,mean=3.082, sd=0.840)),nrow=T,ncol=S)
v.s <- array(rnorm(S*T*K, mean=0,sd=1), dim=c(T,S,K))

## Estimates from Table IV of BLP -- used for comparison and testing
est.blp <- list(alpha=43.501, sigma=c(3.612, 4.628, 1.818, 1.050,
                                      2.056),
                beta=c(-7.061, 2.883, 1.521, -0.122, 3.460),
                gamma=c(0.952, 0.477, 0.619, -.415, -.049, .019))

## Put data into more convenient structure for estimation
est.data <- list(x=as.matrix(cbind(1,BLP[,c("hpwt","air","mpd","space")])),
                 w=as.matrix(cbind(1,log(BLP$hpwt), BLP$air, log(BLP$mpg),
                                   log(BLP$space), BLP$trend)))
est.data$zd <- as.matrix(cbind(est.data$x, Z))
est.data$zs <- as.matrix(cbind(est.data$w, Wiv))
est.data$log.y <- log(y.s)
est.data$log.yp <- list()

## BLP uses log(y - p) in utility function, but some vehicles have
## p>y for some people. BLP do not say what they did in this cases.
## I will take a first order Taylor expansion of log to the left of
## some small number to get an almost log function that is defined
## everywhere. This is very arbitrary though ....
x0 <- 0.1 ## take linear taylor approx to log around x0 for y-p<x0 to
logx0 <- log(x0)
slope <- 1/x0
## avoid log(-)
my.log <- function(x)  ifelse(x>=x0, log(x), logx0 + slope*(x-x0))
dmy.log <- function(x) ifelse(x>=x0, 1/x, slope)
for (t in 1:T) {
  yp <- outer(drop(y.s[t,]), BLP$price[BLP$cdid==t],
              function(x,y) x-y)
  est.data$log.yp[[t]] <- my.log(yp)
  est.data$dlog.yp[[t]] <- dmy.log(yp)
}

## Testing of delta.fn 
t <- 1
inc <- BLP$cdid==t
delta <- rnorm(n=length(BLP$price[inc]))
s <- share.fn(delta, x=drop(est.data$x[inc,]),
              log.y=drop(est.data$log.y[t,]),
              log.yp=est.data$log.yp[[t]],
              v=drop(v.s[t,,]),
              alpha=est.blp$alpha,
              sigma=est.blp$sigma)
d.check <- delta.fn(s, x=drop(est.data$x[inc,]),
                    log.y=drop(est.data$log.y[t,]),
                    log.yp=est.data$log.yp[[t]],
                    v=drop(v.s[t,,]),
                    alpha=est.blp$alpha,
                    sigma=est.blp$sigma, max.iter=1000)
summary((abs(delta-d.check)))

dshare.dp <- function(delta,x,log.y, log.yp, dlog.yp, v,alpha,sigma)
{
  stop("TODO: compute dshare/dp")
  return(dshare)
}

## Testing dshare.dp
if (!require(numDeriv)) install.packages("numDeriv")
library(numDeriv)
t <- 1
dshare.num <- jacobian(function(p) {
  yp <- outer(drop(y.s[t,]), p,
              function(x,y) x-y)
  share.fn(delta, x=drop(est.data$x[inc,]),
           log.y=drop(est.data$log.y[t,]),
           log.yp=my.log(yp),
           v=drop(v.s[t,,]),
           alpha=est.blp$alpha,
           sigma=est.blp$sigma)
}, x=BLP$price[inc])
dshare <- dshare.dp(delta, x=drop(est.data$x[inc,]),
                    log.y=drop(est.data$log.y[t,]),
                    log.yp=est.data$log.yp[[t]], 
                    dlog.yp=est.data$dlog.yp[[t]] ,
                    v=drop(v.s[t,,]),
                    alpha=est.blp$alpha,
                    sigma=est.blp$sigma)
summary(as.vector(dshare-dshare.num))

moments <- function(alpha,sigma, s,p,x,log.yp, log.y,dlog.yp,
                    v,w,zd,zs, W,
                    market.id,firm.id, delta.tol=1e-10,
                    max.iter=100)
{
  ## Find delta and omega for each market
  delta <- rep(NA, length(s))
  omega <- rep(NA, length(s))
  for (t in unique(market.id)) {
    inc <- market.id==t
    delta[inc] <- delta.fn(s[inc], x[inc,],drop(log.y[t,]),
                           log.yp[[t]] ,
                           v=drop(v[t,,]),
                           alpha,sigma, tol=delta.tol,
                           tol.s=1e-6, max.iter=max.iter)
    dShare <- dshare.dp(delta[inc], x[inc,], drop(log.y[t,]),
                        log.yp[[t]], dlog.yp[[t]] , drop(v[t,,]), alpha,sigma)
    dShare <- dShare* outer(firm.id[inc],firm.id[inc],
                            function(x,y) x==y)
    b <- solve(dShare) %*% s[inc]
    omega[inc] <- log(p[inc]-b)
  }
  
  ## Solve for beta and gamma
  X <- as.matrix(rbind(cbind(x,0*w), cbind(0*x,w)))
  Y <- c(delta, omega)
  Z <- as.matrix(rbind(cbind(zd,0*zs), cbind(0*zd, zs)))
  B <- solve(t(X) %*% Z %*% W %*% t(Z) %*% X) %*%
    (t(X) %*% Z %*% W %*% t(Z) %*% Y)
  beta <- B[1:ncol(x)]
  gamma <- B[(1+ncol(x)):(ncol(x)+ncol(w))]
  
  ## Compute GMM objective
  E <- Y - X %*% B
  g <- drop(E)*Z
  G <- colMeans(g)
  obj <- nrow(g)*t(G) %*% W %*% G
  return(list(obj=obj, beta=beta, gamma=gamma))
}

Rprof("blp.prof", line.profiling=TRUE) # start the profiler
q <- moments(est.blp$alpha,est.blp$sigma,s=BLP$share,p=BLP$price,
             x=est.data$x, log.y=est.data$log.y,
             log.yp=est.data$log.yp,
             dlog.yp=est.data$dlog.yp,
             v=v.s,w=est.data$w,zd=est.data$zd,zs=est.data$zs,
             W=diag(1,nrow=(ncol(est.data$zd)+ncol(est.data$zs))),
             market.id=BLP$cdid, firm.id=BLP$firm.id)
Rprof(NULL) # stop the profiler

summaryRprof("blp.prof", lines="both") # show the results