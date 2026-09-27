#!/usr/bin/env Rscript
library(coda)
library(rjags)
library(runjags)
library(tidyverse)
library(loo)

long.time <- read.csv("long.data_new600.csv")
first.tt <- long.time[,2]
last.tt <- long.time[,3]

####time of first visit and last visit#######
N<-length(last.tt)
#participant ID
id<-rep(1:N)
length(id)

t<-round(first.tt)
tt<-round(last.tt)

X1 <- c(rep(1, floor(N/2)), rep(0, N - floor(N/2)))
k.pa<-(tt-t)*4

alpha = c(1,1)

set.seed(123)

#############################################################
X <- as.matrix(read.csv(list.files(pattern="mixJM.X_newdata7.")))
Y <- as.matrix(read.csv(list.files(pattern="mixJM.Y_newdata7.")))
simdat.pe00 <- as.data.frame(read.csv(list.files(pattern="mixJM.rec_newdata7.")))
#############################################################

tt<-tt-0.25
timeS <- as.data.frame(cbind(id,t)) ## left truncation time
timeE <- as.data.frame(cbind(id,tt))

simdat.pe0 <- merge(simdat.pe00, timeS,all=TRUE)
simdat.pe <- subset(simdat.pe0, stop >= t)
simdat.pe <- simdat.pe %>% arrange(id, stop)

N <- length(tt)

# Event times only (status==1)
ev_list <- vector("list", N)

for (i in 1:N) {
  ev_list[[i]] <- simdat.pe$stop[
    simdat.pe$id == i & simdat.pe$status == 1
  ]
}

k.pe <- lengths(ev_list)
max.count <- max(k.pe)

Ti <- t(vapply(ev_list, function(v) {
  if (length(v) == 0) {
    rep(NA_real_, max.count)
  } else {
    c(v, rep(NA_real_, max.count - length(v)))
  }
}, numeric(max.count)))

E <- matrix(0L, nrow=N, ncol=max.count)

for (i in 1:N) {
  if (k.pe[i] > 0)
    E[i, 1:k.pe[i]] <- 1L
}

Ti2 <- Ti
Ti2[is.na(Ti2)] <- 1
#################Readingin data for X, t0, tau vectors#############################
time.t0 <- t
time.tau <- tt

alpha.r = c(1,1)

############Model in the JAGS format#####################
modelrancp <- "
data { 
  for(i in 1:N){
    zeros[i] <- 0
  }
}
model { 
  # -------------------------
  # Subject loop
  # -------------------------
  for(i in 1:N){ 

    # ---- PA (binary) ----
    for(j in 1:k.pa[i]){
      Y[i,j] ~ dbin(p2[i,j], 1)
      p2[i,j] <- p[i,j,z[i]]

      logit(p[i,j,1]) <- c10
        + c[1] * (X[i,j] - cp1[i])
        + c[2] * (X[i,j] - cp1[i]) * (2*step(X[i,j] - cp1[i]) - 1)
        + c[3] * (X[i,j] - cp2[i]) * (2*step(X[i,j] - cp2[i]) - 1)
        + c[4] * X1[i]
        + u1[i]

      logit(p[i,j,2]) <- c20
        + (c[1] - c[2] - c[3]) * X[i,j]
        + c[4] * X1[i]
        + u2[i]
    }

    # Latent class for PA
    z[i] ~ dcat(pi[1:2])

    # Random effects for PA
    u1[i] ~ dnorm(0, u.tau1)
    u2[i] ~ dnorm(0, u.tau2)
    cp1[i] ~ dnorm(cp1.mu, cp1.tau)T(, 21)
    cp2[i] ~ dunif(cp1[i], 21)

    # Center cp1 inside v1
    cp1c[i] <- cp1[i] - cp1.mu
    cp2.mu[i] <- 0.5 * (cp1[i] + 21)
    cp2c[i] <- cp2[i] - cp2.mu[i]
    
    # PA likelihood contribution
    L.a[i]  <- prod( (p2[i,1:k.pa[i]]^Y[i,1:k.pa[i]]) * ((1-p2[i,1:k.pa[i]])^(1-Y[i,1:k.pa[i]])) )
    ll.a[i] <- log(L.a[i])

    # ---- PE (NHPP / Weibull process) ----
    w1[i] ~ dnorm(0, w.tau1)
    w2[i] ~ dnorm(0, w.tau2)
    #v1[i] <- exp(ga10*u1[i] + w1[i] + ga11*cp1c[i] + ga12*cp2c[i]) 
    #v2[i] <- exp(ga20*u2[i] + w2[i])
    
    eta.shared[i] <-
        equals(z[i], 1) *
          (ga10*u1[i] + ga11*cp1c[i] + ga12*cp2c[i]) +
        equals(z[i], 2) *
          (ga20*u2[i])
        
    v1[i] <- exp(w1[i] + eta.shared[i])
    v2[i] <- exp(w2[i] + eta.shared[i])

    # Baseline intensity pieces at event times
    for(j in 1:max.count){
      lambda10[i,j] <- a1 * (Ti2[i,j])^(a1-1)
      lambda20[i,j] <- a2 * (Ti2[i,j])^(a2-1)
      lambda1[i,j] <- lambda10[i,j] * v1[i] * exp(b10 + b[1]*X1[i])
      lambda2[i,j] <- lambda20[i,j] * v2[i] * exp(b20 + b[2]*X1[i])

      loghaz1[i,j] <- E[i,j] * log(lambda1[i,j])
      loghaz2[i,j] <- E[i,j] * log(lambda2[i,j])
    }

    # NHPP log-likelihood for each component:
    # sum log lambda(t_j) - b*(tau^a - t0^a)
    # where b = v*exp(b0+bX)
    logL1[i] <- sum(loghaz1[i,1:max.count]) -
      v1[i] * exp(b10 + b[1]*X1[i]) * (time.tau[i]^a1 - time.t0[i]^a1)

    logL2[i] <- sum(loghaz2[i,1:max.count]) -
      v2[i] * exp(b20 + b[2]*X1[i]) * (time.tau[i]^a2 - time.t0[i]^a2)

    # Mixture over PE components
    z.r[i] ~ dcat(pi.r[1:2]) 
    ll.e[i] <- equals(z.r[i],1) * logL1[i] + equals(z.r[i],2) * logL2[i]

    # zeros trick for custom likelihood
    phi[i] <- max(-ll.e[i] + 10000, 0)
    zeros[i] ~ dpois(phi[i])
  }

  # -------------------------
  # Global summaries
  # -------------------------
  log_lik0.a <- sum(ll.a[])
  log_lik0.e <- sum(ll.e[])
  dev.a <- -2*log_lik0.a
  dev.e <- -2*log_lik0.e

  # -------------------------
  # Priors
  # -------------------------
  pi[1:2]   ~ ddirch(alpha[])
  pi.r[1:2] ~ ddirch(alpha.r[])

  # PA intercept ordering: c20 > c10
  c20_raw ~ dnorm(0, 0.01)
  delta_c ~ dnorm(0, 0.01) T(0,)
  c10 <- c20_raw - delta_c
  c20 <- c20_raw

  for (k in 1:4){
    c[k] ~ dnorm(0,0.01)
  }
  B1 <- c[1] - c[2] - c[3]
  B2 <- c[1] + c[2] - c[3]
  B3 <- c[1] + c[2] + c[3]

  u.tau1 ~ dgamma(16,4) # u.tau1 ~ dgamma(0.001,0.001)
  u.tau.inv1 <- 1/u.tau1

  u.tau2 ~ dgamma(16,4) # u.tau2 ~ dgamma(0.001,0.001)
  u.tau.inv2 <- 1/u.tau2

  cp1.mu ~ dnorm(0, 0.01)
  cp1.tau ~ dgamma(1, 1)
  cp1.tau.inv <- 1/cp1.tau

  a1 ~ dgamma(0.1, 0.1)
  a2 ~ dgamma(0.1, 0.1)
  # PE ordering: b10 > b20
  b20_raw ~ dnorm(0, 0.01)
  delta_b ~ dnorm(0, 0.01) T(0,)
  b10 <- b20_raw + delta_b
  b20 <- b20_raw

  for (p in 1:2){
    b[p] ~ dnorm(0,0.01)
  }

  ga10 ~ dnorm(0, 1)
  ga20 ~ dnorm(0, 1)
  ga11 ~ dnorm(0, 1)
  ga12 ~ dnorm(0, 1)

  w.tau1 ~ dgamma(4, 0.16) # w.tau1 ~ dgamma(2,2)
  w.tau.inv1 <- 1/w.tau1

  w.tau2 ~ dgamma(4, 0.16) # w.tau2 ~ dgamma(2,2)
  w.tau.inv2 <- 1/w.tau2
}
"
####Observed DATA
data <- dump.format(list(N=N, X=X, Y=Y, X1=X1,k.pa=k.pa,max.count=max.count, time.t0=time.t0, time.tau=time.tau, Ti2=Ti2, E=E, alpha=alpha, alpha.r=alpha.r)) 
###initial Values
inits1 <- dump.format(list(
  c20_raw=-2, delta_c=2, c=c(0.4,0.2,0.3,-0.05),
  pi=c(0.5,0.5), pi.r=c(0.5,0.5),
  u.tau1=4, u.tau2=4, cp1.mu=9, cp1.tau=1,
  b20_raw=-2, delta_b=0.5, b=c(0.2,0.3),
  a1=1.1, a2=0.9, w.tau1=25, w.tau2=25,
  ga10=0.3, ga20=0.1, ga11=-0.2, ga12=0.1,
  .RNG.name="base::Super-Duper", .RNG.seed=1
))

inits2 <- dump.format(list(
  c20_raw=-2.1, delta_c=2.1, c=c(0.4,0.2,0.3,-0.05)+0.01,
  pi=c(0.51,0.49), pi.r=c(0.51,0.49),
  u.tau1=3.6, u.tau2=4.4, cp1.mu=9.1, cp1.tau=0.9,
  b20_raw=-2.1, delta_b=0.4, b=c(0.25,0.35),
  a1=1.2, a2=0.85, w.tau1=24, w.tau2=26,
  ga10=0.31, ga20=0.11, ga11=-0.21, ga12=0.11,
  .RNG.name="base::Super-Duper", .RNG.seed=2
))

#### Run the model and produce plots
res <- run.jags(model=modelrancp, adapt = 10000, burnin=10000, sample=6000,  
                monitor=c("B1","B2","B3","c10", "c20","c", "cp1",
                          "pi","pi.r","u.tau.inv1","u.tau.inv2", 
                          "cp1.mu","cp1.tau.inv",
                          "b10","b20","b", "a1","a2","ga10","ga20","ga11","ga12",
                          "w.tau.inv1","w.tau.inv2",
                          "ll.a","ll.e","cp2","cp2.mu","z","z.r"), 
                data=data, n.chains=2, method = "parallel", inits=c(inits1,inits2), thin=10)

summary <- summary(res)
summary
result_df <- as.data.frame(summary)
text <- list.files(pattern="mixJM.X_newdata7.")
num <- unlist(lapply(strsplit(text,'.',fixed=TRUE),function(x) x[[3]]))
write.csv(result_df, paste0("mixJM.newresult11.",num,".csv"))

res_jm <- res$mcmc

## =========================================================
## Helper functions
## =========================================================

colVars <- function(a){
  diff <- a - matrix(colMeans(a), nrow(a), ncol(a), byrow = TRUE)
  colSums(diff^2) / (nrow(a) - 1)
}

log_mean_exp <- function(x){
  m <- max(x)
  m + log(mean(exp(x - m)))
}

waic_from_loglik <- function(log_lik){
  
  lppd_i <- apply(log_lik, 2, log_mean_exp)
  lppd <- sum(lppd_i)
  
  p_waic_1 <- 2 * sum(lppd_i - colMeans(log_lik))
  p_waic_2 <- sum(colVars(log_lik))
  
  WAIC <- -2 * (lppd - p_waic_2)
  
  list(
    lppd = lppd,
    p_waic_1 = p_waic_1,
    p_waic_2 = p_waic_2,
    WAIC = WAIC
  )
}

dic_from_loglik <- function(log_lik){
  
  D <- -2 * rowSums(log_lik)
  
  mean_deviance <- mean(D)
  pD <- var(D) / 2
  DIC <- mean_deviance + pD
  
  list(
    mean_deviance = mean_deviance,
    pD = pD,
    DIC = DIC
  )
}


## =========================================================
## Combine posterior draws across chains
## =========================================================

post <- do.call(rbind, lapply(res_jm, as.matrix))
cn <- colnames(post)


## =========================================================
## Locate and order subject-level log-likelihood columns
## =========================================================

idx.a <- grep("^ll\\.a\\[[0-9]+\\]$", cn)
idx.e <- grep("^ll\\.e\\[[0-9]+\\]$", cn)

if(length(idx.a) == 0)
  stop("No ll.a[i] columns found in posterior samples.")

if(length(idx.e) == 0)
  stop("No ll.e[i] columns found in posterior samples.")

idx.a <- idx.a[
  order(as.integer(
    sub("^ll\\.a\\[([0-9]+)\\]$", "\\1", cn[idx.a])
  ))
]

idx.e <- idx.e[
  order(as.integer(
    sub("^ll\\.e\\[([0-9]+)\\]$", "\\1", cn[idx.e])
  ))
]

if(length(idx.a) != length(idx.e))
  stop("Different number of ll.a[i] and ll.e[i] columns.")


## =========================================================
## Subject-level log-likelihood matrices
## =========================================================

ll.a.mat <- post[, idx.a, drop = FALSE]
ll.e.mat <- post[, idx.e, drop = FALSE]

ll.total.mat <- ll.a.mat + ll.e.mat


## =========================================================
## DIC
## =========================================================

dic.pa    <- dic_from_loglik(ll.a.mat)
dic.pe    <- dic_from_loglik(ll.e.mat)
dic.total <- dic_from_loglik(ll.total.mat)


## =========================================================
## WAIC
## =========================================================

waic.pa    <- waic_from_loglik(ll.a.mat)
waic.pe    <- waic_from_loglik(ll.e.mat)
waic.total <- waic_from_loglik(ll.total.mat)


## =========================================================
## PSIS-LOO
## ========================================================= 

loo.pa    <- loo(ll.a.mat)
loo.pe    <- loo(ll.e.mat)
loo.total <- loo(ll.total.mat)


## =========================================================
## Final output: DIC + WAIC + LOO only
## =========================================================

dicwaic_df <- data.frame(
  DIC_PA    = dic.pa$DIC,
  WAIC_PA   = waic.pa$WAIC,
  LOOIC_PA  = -2 * loo.pa$estimates["elpd_loo", "Estimate"],
  
  DIC_PE    = dic.pe$DIC,
  WAIC_PE   = waic.pe$WAIC,
  LOOIC_PE  = -2 * loo.pe$estimates["elpd_loo", "Estimate"],
  
  DIC_Total   = dic.total$DIC,
  WAIC_Total  = waic.total$WAIC,
  LOOIC_Total = -2 * loo.total$estimates["elpd_loo", "Estimate"]
)

write.csv(
  dicwaic_df,
  paste0("dicwaic_mixJM11.", num, ".csv"),
  row.names = FALSE
)


post <- do.call(rbind, lapply(res_jm, as.matrix))
cn <- colnames(post)

get_subject_id <- function(x, prefix){
  as.integer(sub(paste0("^", prefix, "\\[([0-9]+)\\]$"), "\\1", x))
}

idx.z   <- grep("^z\\[[0-9]+\\]$", cn)
idx.cp1 <- grep("^cp1\\[[0-9]+\\]$", cn)
idx.cp2 <- grep("^cp2\\[[0-9]+\\]$", cn)

idx.z   <- idx.z[order(get_subject_id(cn[idx.z], "z"))]
idx.cp1 <- idx.cp1[order(get_subject_id(cn[idx.cp1], "cp1"))]
idx.cp2 <- idx.cp2[order(get_subject_id(cn[idx.cp2], "cp2"))]

if(length(idx.z) != N) stop("Number of z[i] parameters does not equal N.")
if(length(idx.cp1) != N) stop("Number of cp1[i] parameters does not equal N.")
if(length(idx.cp2) != N) stop("Number of cp2[i] parameters does not equal N.")

z.mat   <- post[, idx.z, drop=FALSE]
cp1.mat <- post[, idx.cp1, drop=FALSE]
cp2.mat <- post[, idx.cp2, drop=FALSE]

## Classification: Component 1 if P(z=1)>=0.75, otherwise Component 2
prob.comp1 <- colMeans(z.mat == 1)
component <- ifelse(prob.comp1 >= 0.75, 1, 2)

## Subject-specific posterior means
cp1.post.mean <- colMeans(cp1.mat)
cp2.post.mean <- colMeans(cp2.mat)

classify.dat <- data.frame(
  ID=1:N,
  Prob_Component1=prob.comp1,
  Component=component,
  cp1_post_mean=cp1.post.mean,
  cp2_post_mean=cp2.post.mean
)

## Number and percentage in each component
n.comp1 <- sum(component == 1)
n.comp2 <- sum(component == 2)
pct.comp1 <- 100*n.comp1/N
pct.comp2 <- 100*n.comp2/N

component1.id <- which(component == 1)
component2.id <- which(component == 2)

## cp1 and cp2 summaries among Component 1 subjects
if(length(component1.id) > 0){
  cp1.subject.mean <- cp1.post.mean[component1.id]
  cp2.subject.mean <- cp2.post.mean[component1.id]
  
  cp1.comp1.mean <- mean(cp1.subject.mean)
  cp1.comp1.lower <- as.numeric(quantile(cp1.subject.mean, 0.025, na.rm=TRUE))
  cp1.comp1.upper <- as.numeric(quantile(cp1.subject.mean, 0.975, na.rm=TRUE))
  
  cp2.comp1.mean <- mean(cp2.subject.mean)
  cp2.comp1.lower <- as.numeric(quantile(cp2.subject.mean, 0.025, na.rm=TRUE))
  cp2.comp1.upper <- as.numeric(quantile(cp2.subject.mean, 0.975, na.rm=TRUE))
} else {
  cp1.comp1.mean <- cp1.comp1.lower <- cp1.comp1.upper <- NA_real_
  cp2.comp1.mean <- cp2.comp1.lower <- cp2.comp1.upper <- NA_real_
}

classification_summary <- data.frame(
  N=N,
  N_Component1=n.comp1,
  Percent_Component1=pct.comp1,
  N_Component2=n.comp2,
  Percent_Component2=pct.comp2,
  Mean_cp1_Component1=cp1.comp1.mean,
  Lower95_cp1_Component1=cp1.comp1.lower,
  Upper95_cp1_Component1=cp1.comp1.upper,
  Mean_cp2_Component1=cp2.comp1.mean,
  Lower95_cp2_Component1=cp2.comp1.lower,
  Upper95_cp2_Component1=cp2.comp1.upper
)


write.csv(
  classification_summary,
  paste0("mixJM.classification_summary11.", num, ".csv"),
  row.names=FALSE
)