#!/usr/bin/env Rscript
library(coda)
library(rjags)
library(runjags)
library(tidyverse)
library(mcmcplots)

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
X <- as.matrix(read.csv(list.files(pattern="mixJM.X_newdata1.")))
Y <- as.matrix(read.csv(list.files(pattern="mixJM.Y_newdata1.")))
simdat.pe00 <- as.data.frame(read.csv(list.files(pattern="mixJM.rec_newdata1.")))
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
        + c[3] * X1[i]
        + u1[i]

      logit(p[i,j,2]) <- c20
        + (c[1] - c[2]) * X[i,j]
        + c[3] * X1[i]
        + u2[i]
    }

    # Latent class for PA
    z[i] ~ dcat(pi[1:2])

    # Random effects for PA
    u1[i] ~ dnorm(0, u.tau1)
    u2[i] ~ dnorm(0, u.tau2)
    cp1[i] ~ dnorm(cp1.mu, cp1.tau)

    # Center cp1 inside v1
    cp1c[i] <- cp1[i] - cp1.mu
    
    # PA likelihood contribution
    L.a[i]  <- prod( (p2[i,1:k.pa[i]]^Y[i,1:k.pa[i]]) * ((1-p2[i,1:k.pa[i]])^(1-Y[i,1:k.pa[i]])) )
    ll.a[i] <- log(L.a[i])

    # ---- PE (NHPP / Weibull process) ----
    w1[i] ~ dnorm(0, w.tau1)
    w2[i] ~ dnorm(0, w.tau2)
    v1[i] <- exp(ga10*u1[i] + w1[i] + ga11*cp1c[i]) # v1[i] <- exp(ga10*u1[i] + w1[i] + ga11*cp1[i])
    v2[i] <- exp(ga20*u2[i] + w2[i])

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
    phi[i] <- max(-ll.e[i] + 1000000, 0)
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

  for (k in 1:3){
    c[k] ~ dnorm(0,0.01)
  }
  B1 <- c[1] - c[2]
  B2 <- c[1] + c[2]

  u.tau1 ~ dgamma(16,4) # u.tau1 ~ dgamma(0.001,0.001)
  u.tau.inv1 <- 1/u.tau1

  u.tau2 ~ dgamma(16,4) # u.tau2 ~ dgamma(0.001,0.001)
  u.tau.inv2 <- 1/u.tau2

  cp1.mu ~ dnorm(0,0.01)
  cp1.tau ~ dgamma(1,1)
  cp1.tau.inv <- 1/cp1.tau

  a1 ~ dgamma(0.01,0.01)
  a2 ~ dgamma(0.01,0.01)
  # PE ordering: b10 < b20
  b20_raw ~ dnorm(0, 0.01)
  delta_b ~ dnorm(0, 0.01) T(0,)
  b10 <- b20_raw - delta_b
  b20 <- b20_raw

  for (p in 1:2){
    b[p] ~ dnorm(0,0.01)
  }

  ga10 ~ dnorm(0,1)
  ga20 ~ dnorm(0,1)
  ga11 ~ dnorm(0,1)

  w.tau1 ~ dgamma(25, 2.25) # w.tau1 ~ dgamma(2,2)
  w.tau.inv1 <- 1/w.tau1

  w.tau2 ~ dgamma(25, 2.25) # w.tau2 ~ dgamma(2,2)
  w.tau.inv2 <- 1/w.tau2
}
"

####Observed DATA
data <- dump.format(list(N=N, X=X, Y=Y, X1=X1,k.pa=k.pa,max.count=max.count, time.t0=time.t0, time.tau=time.tau, Ti2=Ti2, E=E, alpha=alpha, alpha.r=alpha.r)) 
###initial Values
inits1 <- dump.format(list(c20=-2.6, delta_c=0.7, c=c(0.3,0.3,-0.05), pi=c(0.55,0.45), pi.r=c(0.6,0.4), u.tau1=4,u.tau2=4, cp1.mu=14, cp1.tau=1,
                           b20_raw=-2, delta_b=2, b=c(0.2,0.3), a1=1.8,a2=1.6, w.tau1=11.1, w.tau2=11.1, ga10=1, ga20=-0.2, ga11=-0.1,
                           .RNG.name="base::Super-Duper", .RNG.seed=1)) 
inits2 <- dump.format(list(c20=-2.5, delta_c=0.6, c=c(0.3,0.3,-0.05)+0.01, pi=c(0.56,0.44), pi.r=c(0.59,0.41), u.tau1=3.6,u.tau2=4.4, cp1.mu=14.1, cp1.tau=0.9,
                           b20_raw=-1.9, delta_b=2.2, b=c(0.2,0.3)+0.1, a1=1.75,a2=1.55, w.tau1=10, w.tau2=12, ga10=1.1, ga20=-0.1, ga11=-0.08,
                           .RNG.name="base::Super-Duper", .RNG.seed=2))

#### Run the model and produce plots
res <- run.jags(model=modelrancp, burnin=10000, sample=5000,  
                monitor=c("B1","B2","c10", "c20","c", "cp1",
                          "pi","pi.r","u.tau.inv1","u.tau.inv2", "u.tau1","u.tau2",
                          "cp1.mu","cp1.tau.inv","cp1.tau",
                          "b10","b20","b", "a1","a2","ga10","ga20","ga11",
                          "w.tau1","w.tau2","w.tau.inv1","w.tau.inv2","c20_raw", "delta_c","b20_raw","delta_b",
                          "ll.a","ll.e","dev.a","dev.e"), 
                data=data, n.chains=2, method = "parallel", inits=c(inits1,inits2), thin=15)

summary <- summary(res)
summary
result_df <- as.data.frame(summary)
text <- list.files(pattern="mixJM.X_newdata1.")
num <- unlist(lapply(strsplit(text,'.',fixed=TRUE),function(x) x[[3]]))
write.csv(result_df, paste0("mixJM.newresult5.",num,".csv"))
save(res, file=paste0("mixJM.newres5.",num,".RData"))

res_jm <- res$mcmc
vars<-mcmc.list(res_jm[[1]][,c(1:16)],res_jm[[2]][,c(1:16)])
pdf(file = paste0("mixJM.newtraceplot5.",num,".pdf"),   # The directory you want to save the file in
    width = 4, # The width of the plot in inches
    height = 4) # The height of the plot in inches
traplot(vars)
dev.off()
