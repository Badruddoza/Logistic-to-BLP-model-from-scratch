# BLP using cereal data
## following Tyler Mangin
while(!require(AER)){install.packages("AER")}
while(!require(SQUAREM)){install.packages("SQUAREM")}
while(!require(BB)){install.packages("BB")}
while(!require(plyr)){install.packages("plyr")}

share.fld =     "share"
prod.id.fld =   "brand"
mkt.id.fld =    "city"
prc.fld =       "price"
x.var.flds =    c("sugar",
                  "mushy")

#Set up data
dat <- dat[dat[, share.fld] > 0, ]
dat <- dat[order(dat[, mkt.id.fld], dat[, prod.id.fld]), ]
JM <- nrow(dat)

#Number of characteristics (including constant and price)
X <- as.matrix(cbind(ones = rep(1, JM), dat[, c(x.var.flds, prc.fld)]));
K <- ncol(X)

#market object
mkt.id <- dat[, mkt.id.fld];
#shares object
s.jm <- as.vector(dat[, share.fld]);
#sum of product shares by market
temp <- aggregate(s.jm, by = list(mkt.id = mkt.id), sum);

summary(temp)

#Compute the outside good market share by market
sum1 <- temp$x[match(mkt.id, temp$mkt.id)];
s.j0 <- as.vector(1 - sum1);
rm(temp, sum1);
### LOGIT MODEL
# delta object
dat[, "delta"] <- Y <- log(s.jm) - log(s.j0);

#OLS
fm.olsreg = paste0("delta ~ ", 
                   paste(x.var.flds, collapse = " + "), " + ", 
                   paste(prc.fld, collapse = " + "))
ols = lm(data = dat,
         formula = fm.olsreg)
summary(ols)

#2SLS
prc.iv.flds = c("z1",
                "z2",
                "z3")
beta.est = NULL;

str.ivreg.y <- "delta ~ "
str.ivreg.x <- paste(x.var.flds, collapse = " + ")
str.ivreg.prc <- paste(prc.fld, collapse = " + ")
str.ivreg.iv <- paste(prc.iv.flds, collapse = " + ")
print("2SLS specification:")
print(fm.ivreg <- paste0(str.ivreg.y, str.ivreg.x, " + ", str.ivreg.prc, " | ", str.ivreg.x, " + ", str.ivreg.iv))
rm(str.ivreg.y, str.ivreg.x, str.ivreg.prc, str.ivreg.iv)

print("2SLS beta estimate:")
print(summary(mo.ivreg <- ivreg(fm.ivreg, data = dat, x = TRUE)))
beta.est <- summary(mo.ivreg)$coef[, 1:2]
#Z = instrumental variable matrix include exogenous X's
Z <- as.matrix(mo.ivreg$x$instruments)
PZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z);
theta1 <- coef(mo.ivreg);
xi.hat <- as.vector(mo.ivreg$resid);
Z.hat <- Z * matrix(rep(xi.hat, ncol(Z)), ncol = ncol(Z))
W.inv <- try(solve(t(Z.hat) %*% Z.hat), silent = FALSE)
if("matrix" == class(W.inv)){
  PZ <- Z %*% W.inv %*% t(Z);
  PX.inv <- solve(t(X) %*% PZ %*% X)
  theta1 <- PX.inv %*% t(X) %*% PZ %*% Y
  xi.hat <- Y - X %*% theta1
  X.hat <- (PZ %*% X) * matrix(rep(xi.hat, K), ncol = K)
  tsls.se <- sqrt(diag(PX.inv %*% t(X.hat) %*% X.hat %*% PX.inv))
  # print("GMM step 2 updated theta1 estimate:")
  # print(beta.est <- data.frame(beta.est = theta1, se.est = tsls.se))
}
dat[, "xi.hat"] <- xi.hat

##Random Coefficents
## Matrix of individuals' characteristics ##
#number of simulated consumers
n.sim = 100
#Standard normal distribution draws, one for each characteristic in X
#columns are simulated consumers, rows are variables in X (including constant and price)
v = matrix(rnorm(K * n.sim), nrow = K, ncol = n.sim)
##Optimization routine
multiStart(par,
           fn,
           gr,
           action = c("solve", "optimize"),
           method=c(2,3,1),
           lower=-Inf,
           upper=Inf,
           project=NULL,
           projectArgs=NULL,
           control=list(),
           quiet=FALSE,
           details=FALSE)
ind_sh <- function(delta.in, mu.in){
  # This function computes the "individual" probabilities of choosing each brand
  # Requires global variables: mkt.id, X, v
  numer <- exp(mu.in) * matrix(rep(exp(delta.in), n.sim), ncol = n.sim);
  denom <- as.matrix(do.call("rbind", lapply(mkt.id, function(tt){
    1 + colSums(numer[mkt.id %in% tt, ])
  })))
  return(numer / denom);  
}
blp_inner <- function(delta.in, mu.in) {
  # Computes a single update of the BLP (1995) contraction mapping.
  # of market level predicted shares.
  # This single-update function is required by SQUAREM, see Varadhan and
  # Roland (SJS, 2008), and Roland and Varadhan (ANM, 2005)
  # INPUT
  #   delta.in : current value of delta vector
  #   mu.in: current mu matrix
  # Requires global variables: s.jm
  # OUTPUT
  #   delta.out : delta vector that equates observed with predicted market shares
  pred.s <- rowMeans(ind_sh(delta.in, mu.in));
  delta.out <- delta.in + log(s.jm) - log(pred.s)
  return(delta.out)
}

gmm_obj <- function(theta2){
  # This function computes the GMM objective function
  # Requires global variable inputs: X, v, delta, a, W
  # Outputs: theta1, xi.hat
  print(paste0("GMM Loop number: ", Sys.time()))
  print(a <<- a + 1)
  print("Updated theta2 estimate:")
  print(theta2)
  print("Change in theta2 estimate:")
  print(theta.chg <- as.numeric(theta2 - theta2.prev));
  if(sum(theta.chg != 0) <= 2){
    delta <- dat[, "delta"];
  } else {
    delta <- Y;
  }
  theta2.prev <<- theta2;
  
  mu <- X %*% diag(theta2) %*% v;
  
  print("Running SQUAREM contraction mapping")
  print(system.time(
    squarem.output <- squarem(par = delta, fixptfn = blp_inner, mu.in = mu, control = list(trace = TRUE))
  ));
  delta <- squarem.output$par
  print(summary(dat[, "delta"] - delta));
  dat[, "delta"] <<- delta;
  
  mo.ivreg <- ivreg(fm.ivreg, data = dat, x = TRUE)
  theta1 <<- coef(mo.ivreg);
  xi.hat <<- as.vector(mo.ivreg$resid);
  Z.hat <- Z * matrix(rep(xi.hat, ncol(Z)), ncol = ncol(Z))
  W.inv <- try(solve(t(Z.hat) %*% Z.hat), silent = FALSE)
  
  if("matrix" == class(W.inv)){
    
    PX.inv <- solve(t(X) %*% PZ %*% X)
    theta1 <<- PX.inv %*% t(X) %*% PZ %*% delta
    xi.hat <<- delta - X %*% theta1
    X.hat <- (PZ %*% X) * matrix(rep(xi.hat, K), ncol = K)
    tsls.se <- sqrt(diag(PX.inv %*% t(X.hat) %*% X.hat %*% PX.inv))
    print("GMM step 2 updated theta1 estimate:")
    print(beta.est <<- data.frame(beta.est = theta1, beta.se = tsls.se, sigma.est = theta2))
    print("made it here")
  }
  
  dat[, "xi.hat"] <<- xi.hat
  f <- t(xi.hat) %*% Z %*% W.inv %*% t(Z) %*% xi.hat;

  print("Updated GMM objective:")
  print(f <- as.numeric(f));
  return(f)
}

jacobian <- function(delta.in, theta.in){
  print(paste0("Calculating Jacobian matrix, ", Sys.time()))
  #Requires global variables X, v, mkt.id
  mu1 <- X %*% diag(theta.in) %*% v;
  ind.shares <- ind_sh(delta.in, mu1);
  K <- ncol(X);
  print(paste0("Calculating dsigma matrix, ", Sys.time()))
  dsigma <- lapply(l.Xv, function(x){
    temp2 <- x * ind.shares;
    temp3 <- as.matrix(do.call("rbind", lapply(mkt.id, function(m){
      colSums(temp2[mkt.id %in% m, ])
    })));
    dsigma.res <- rowMeans(temp2 - ind.shares * temp3);
    return(dsigma.res)
  })
  dsigma <- as.matrix(do.call("cbind", dsigma))
  print(paste0("Calculating ddelta matrices, ", Sys.time()))
  ddelta <- list()
  for(m in mkt.id){
    if(m %in% names(ddelta)){next}
    temp1 <- as.matrix(ind.shares[mkt.id %in% m, ]);
    H1 <- temp1 %*% t(temp1);
    H2 <- diag(rowSums(temp1));
    H <- (H2 - H1) / n.sim;
    H.inv <- solve(H);
    ddelta[[as.character(m)]] <- H.inv %*% dsigma[mkt.id %in% m, ];
    rm(temp1, H1, H2, H, H.inv)
  }
  ddelta <- as.matrix(do.call("rbind", ddelta));
  return(ddelta)
}

gradient_obj <- function(theta2){
  #Requires global variables PZ, delta, xi.hat
  print(system.time(jacobian_res <<- jacobian(as.vector(dat[, "delta"]), theta2)))
  print(paste0("Updated gradient:", Sys.time()))
  print(f <- -2 * as.numeric(t(jacobian_res) %*% PZ %*% xi.hat));
  #######
  L <- ncol(Z)
  covg <- matrix(0, nrow = L, ncol = L)
  for(i in 1:JM){
    covg <- covg + (Z[i, ] %*% t(Z[i, ])) * xi.hat[i]^2
  }
  d.delta <- jacobian_res;
  Dg <- t(d.delta) %*% Z
  p.Dg <- try(solve(Dg %*% W.inv %*% t(Dg)))
  cov.mat <- p.Dg %*% (Dg %*% W.inv %*% covg %*% W.inv %*% t(Dg)) %*% p.Dg
  beta.est$sigma.se <<- sqrt(diag(cov.mat));
  print(paste0("Updated coefficients table:", Sys.time()))
  print(beta.est)
  write.csv(beta.est, file = paste0("BLP_beta_est_", Sys.Date(), ".csv"))
  #######
  return(as.numeric(f))
}

#Starting point
print("Sigma guess:")
## [1] "Sigma guess:"
# tsls.se = beta.est[, 2]
print(theta2 <- 0.5 * tsls.se);
##        ones       sugar       mushy       price 
## 0.727581577 0.008837817 0.049185287 6.303916821
theta2 = t(theta2)

theta2.prev <- theta2;

# Break X and v matrices into list variables 
# in attempt to expedite calculation of the Jacobian matrix
l.X <- lapply(1:K, function(k){
  return(X[, k])
})
l.v <- lapply(1:K, function(k){
  return(v[k, ])
})
l.Xv <- lapply(1:K, function(k){
  l.X[[k]] %*% t(l.v[[k]]);
})

print("Estimating random coefficients multinomial logit")
a <- 0;
beta.est <- NULL;
print(system.time(
  theta.est <- multiStart(par = theta2, fn = gmm_obj, gr = gradient_obj, lower = 0, control = list(trace = TRUE), action = "optimize")
));
save(theta.est, file = paste0("theta_est_.RData"))
